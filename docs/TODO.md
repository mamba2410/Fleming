# TODO

## Known bugs/issues

- Commented out line in `FluxFinder.py:262` in `make_light_curves()`. Seems median is zero or negative
- Breaks when finding shifts with l136_0 S4I030, not enough values for Gaussian

Runtime warnings/invalid values on l137_0

- DA:228
- FF:432
- DA:604

## Improvements

- New way of determining variables?
- Take field with obvious variable, implement simple variable detections to cut down dataset
- Need sigma clipping
- Can FluxFinder give errors?
- Complicated way of finding variables:
	- smooth the curve with running optimal average
	- sigma clip
	- smooth again
	- find max and min values
	- reject if opt avg fit has too large standard deviation?
	- if variables overlap with stars that hit edge of field, ignore
- raw/reduced images to be stored in separate directory
- Stack images and use that as the catalogue image?
- Assumes catalogue is first image (where exactly?)
- Thumbnail size in `Constants.py`
- Make a `Catalogue` class which has numpy array of sources to save file reading.
- Docs on the readme about the 'public' functions that users can run.
- Periodogram 

## Possible bugs

- Fix `WARNING: Fit may be unsuccessful` in `ShiftFinder.py`

## Future plans

- Have a log level system to unclutter output.

