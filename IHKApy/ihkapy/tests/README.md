# Tests

This is a procedural codebase and the proceedures a split in four steps that can be summarized as follows:
1. Turn the EDF raw data files into binary wavelet banks
2. Compute feature sets based on recording annotation metadata
3. Train a model to classify the data


Featurising our data runs relatively fast but requires large amounds of data. (40GB per 24h session). So our tests integration tests rely on having a local copy of the data. 
Turning EDFs into wavelet banks is straightforward enough, and only requires small amounts of data, so we can test on much smaller and more portable files. 

