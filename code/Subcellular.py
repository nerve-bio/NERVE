"""Run PSORTb localization and require a result for every input protein."""

import json
import logging
import math
import urllib.error
import urllib.request
from operator import attrgetter

from code.Utils import build_id_index, match_protein


class Localization:
    def __init__(self, localization, reliability):
        self.localization = str(localization)
        self.reliability = float(reliability)


def parse_psortb_result(result, proteins):
    """Parse the existing PSORTb text format without accepting partial output."""
    index = build_id_index(proteins)
    found = {}
    separator = '-------------------------------------------------------------------------------\n\n'
    for entry in result.split(separator):
        if 'SeqID: ' not in entry:
            continue
        identifier = entry.split('SeqID: ', 1)[1].splitlines()[0].strip()
        protein = match_protein(identifier, index, proteins)
        if protein is None:
            raise ValueError(f'PSORTb returned an unknown protein: {identifier}')
        if protein.accession in found:
            raise ValueError(f'PSORTb returned duplicate protein: {identifier}')
        if 'Localization Scores:' not in entry or 'Final' not in entry:
            raise ValueError(f'PSORTb omitted localization scores for {identifier}')
        scores = entry.split('Localization Scores:', 1)[1].split('Final', 1)[0]
        locations = []
        for line in scores.strip().splitlines():
            parts = line.split()
            if len(parts) != 2:
                raise ValueError(f'PSORTb returned a malformed score for {identifier}: {line!r}')
            location = Localization(parts[0], parts[1])
            if not math.isfinite(location.reliability):
                raise ValueError(f'PSORTb returned a nonfinite score for {identifier}')
            locations.append(location)
        if not locations:
            raise ValueError(f'PSORTb returned no localization for {identifier}')
        locations.sort(key=attrgetter('reliability'), reverse=True)
        if locations[0].reliability <= 3.0:
            locations = [Localization('Unknown', 0.0)]
        found[protein.accession] = locations
    missing = [protein.accession for protein in proteins if protein.accession not in found]
    if missing:
        raise ValueError(f'PSORTb omitted {len(missing)} proteins; first: {missing[0]}')
    return found


def psortb(list_of_proteins, working_dir, gram, proteome1) -> list:
    """Submit one proteome to the PSORTb service and attribute every result."""
    with open(proteome1) as handle:
        sequence_text = handle.read()
    payload = json.dumps({'gram': gram, 'seq': sequence_text}).encode('utf-8')
    request = urllib.request.Request('http://psortb:8080/execute', data=payload,
                                     headers={'Content-Type': 'application/json; charset=utf-8'})
    try:
        # A whole large-proteome request can legitimately take over 25 minutes.
        with urllib.request.urlopen(request, timeout=3600) as response:
            data = json.load(response)
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as error:
        raise RuntimeError(f'PSORTb request failed: {error}') from error
    if not isinstance(data, dict) or not isinstance(data.get('result'), str):
        raise RuntimeError('PSORTb returned no usable result text')
    if data.get('error'):
        raise RuntimeError(f'PSORTb service error: {data["error"]}')
    logging.debug('PSORTb response received for %d proteins', len(list_of_proteins))
    found = parse_psortb_result(data['result'], list_of_proteins)
    for protein in list_of_proteins:
        protein.localization = found[protein.accession]
    return list_of_proteins
