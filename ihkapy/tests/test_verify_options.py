"""Verify that the Options.toml is logically consistent with itsself."""

import toml

ops_path = "../Options.toml"
with open(ops_path,"r") as f:
    ops_string = f.read()
ops = toml.loads(ops_string)

### params.data options
data = ops["params"]["data"]
# Make sure NUM_FREQ, etc. are compatible
assert data["N_CHAN_BINARY"] == data["NUM_FREQ"] * len(data["TS_FEATURES"]) + 1, "N_CHAN_BINARY = {data['N_CHAN_BINARY']} must be compatible with num_freq = {data['NUM_FREQ']} and ts_features = {data['TS_FEATURES']}"
# Check that the relative Amplitude and Phase indices are compatible  
# with NUM_FREQ
for fidx in data["AMP_IDX"]:
    assert fidx == int(fidx) and fidx >= 0 and fidx < data["NUM_FREQ"]
for fidx in data["PH_IDX"]:
    assert fidx == int(fidx) and fidx >= 0 and fidx < data["NUM_FREQ"]
NCB = data["N_CHAN_BINARY"]
for fidx in data["AMP_FREQ_IDX_ALL"]:
    assert fidx not in data["PH_FREQ_IDX_ALL"] 
    assert fidx < NCB and fidx > 0
for fidx in data["PH_FREQ_IDX_ALL"]:
    assert fidx not in data["AMP_FREQ_IDX_ALL"]
    assert fidx < NCB and fidx > 0
print("Data tests passed")

### params.feature options
feat = ops["params"]["feature"]
assert feat["N_PREICTAL_BINS"] == len(feat["PREICTAL"]["BINS"]) , "N_PREICTAL_BINS incompatible with len(preictal_bins)"
for i in feat["PREICTAL"]["PCT"]: assert i <= 1.0 and i >= 0.0
assert len(feat["BIN_NAMES"]) == feat["N_PREICTAL_BINS"] + 2
assert feat["POSTICTAL"]["DELAY"] > 0, "postictal delay must be positive"
print("feature tests passed")

### params.model.training options
train = ops["params"]["model"]["training"]
tvt = train["train_val_test"]
assert sum(tvt) == 1.0
print("model training params tests passed")

# TODO: go through entire codebase and find where options are used
# then add an assert here to make sure the thing used will exist
# basically verify all the keys are being used correctly

# It also might be a though to programatically run through codebase
# to verify that the orthodox convention for naming dicts is respected
# and the keys used are all orthodox too. This is probably overkill.



