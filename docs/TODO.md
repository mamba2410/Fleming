# TODO

## Known bugs/issues


## Improvements

- New way of determining variables. Split into `VariableFinder.py`
	- Already have standard deviation search, does not detect low amplitude variables.
	- Implement period search for every light curve.
- Need sigma clipping
- Stack images and use that as the catalogue image?
- Assumes catalogue is first image (where exactly? Something related to shifts)
- Thumbnail size in `Constants.py`
- Make a `Catalogue` class which has numpy array of sources to save file reading.
- Docs on the readme about the 'public' functions that users can run.


## Possible bugs

- Fix `WARNING: Fit may be unsuccessful` in `ShiftFinder.py`
- Runtime warnings/invalid values on l137_0
- Commented out if statement in `FluxFinder.py:222` in `make_light_curves()`.
	Seems median is sometimes zero or negative.

## Future plans

- Have a log level system to unclutter output.

