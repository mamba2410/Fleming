# TODO

## Immediate improvements

- Reformat `Constants.py` -> `Config.py` and have it be an object to pass around
- Convert a lot of arrays to numpy arrays.
- Merge/cache functions and file names
- Thumbnail size in `Constants.py`
- Do something with the WCS (ie add RA/DEC to catalogue)
- Make a `Catalogue` class which has numpy array of sources to save file reading.
- Assumes catalogue is first image.

## Possible bugs

- Fix `WARNING: Fit may be unsuccessful` in `ShiftFinder.py`
- Rerunning something (either in DA or FF) adds more times to the times file.


## Future plans

- Have a log level system to unclutter output.

