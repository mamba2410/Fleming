# TODO

## Immediate improvements

- Reformat `Constants.py` -> `Config.py` and have it be an object to pass around
- Clean up `Constants.py` and make it be used consistently.
	- `_dir` is a full directory path, `_subdir` is a sub-path
	- `_path` is a full path to a file, `_fname` is a base file name
- Make directories more consistent
- Merge/cache functions and file names
- Convert a lot of arrays to numpy arrays.
- Have a log level system to unclutter output.
- Image file path made in `Utilities.py`
- Thumbnail size in `Constants.py`


## Possible bugs

- Fix `WARNING: Fit may be unsuccessful` in `ShiftFinder.py`
- Rerunning something (either in DA or FF) adds more times to the times file.


## Future plans


