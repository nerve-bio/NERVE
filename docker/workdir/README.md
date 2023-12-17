# Test and tuning of NERVE2

Test and tuning was performed on this dataset: [test_antigens_summary_v2.xlsx](../../database/antigens/test_antigens_summary_v2.xlsx). The set was split into test and tuning considering the similarity of the antigens with those used to train other reverse vaccinology tools in order to have an independent dataset for comparing performancies with these tools.
The proteome associated with each antigen was submitted to NERVE2 and analysed with the select module deactivated.
The tuning set was then used to test different variations of the select modules in order to find the best decision tree combination and the best set of parameters.
The following versions of the select module were then tested on the test dataset:

Versions with the final version of select:
- v3b: localization + adhesin + immuno human
- V3i: localization + adhesin + immuno human + razor
- V3l: localization + adhesin + immuno human + razor + mouse
- V3M: localization + adhesin + immuno human + razor + mouse + virulent

Test and tuning datasets for both NERVE1 and NERVE2 can be found in [test](test/) and [tuning](tuning/)