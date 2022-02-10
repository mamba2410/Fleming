# TODO

## Known bugs/issues

- Are adjusted light curves okay? Some seem wild
- Commented out line in `FluxFinder.py:262` in `make_light_curves()`. Seems median is zero or negative

## Improvements

- New way of determining variables?
	- smooth the curve with running optimal average
	- sigma clip
	- smooth again
	- find max and min values
	- reject if opt avg fit has too large standard deviation?
	- if variables overlap with stars that hit edge of field, ignore
- raw/reduced images to be stored in separate directory
- sigma clip light curves?
- concept of categorising light curves? messy (wide std), flat (very unlikely variable), human (unsure)
- Thumbnail size in `Constants.py`
- Do something with the WCS (ie add RA/DEC to catalogue)
- Make a `Catalogue` class which has numpy array of sources to save file reading.
- Assumes catalogue is first image (where exactly?)

## Possible bugs

- Fix `WARNING: Fit may be unsuccessful` in `ShiftFinder.py`

## Future plans

- Have a log level system to unclutter output.

