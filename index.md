<p align="center">
  <img src="https://github.com/user-attachments/assets/478422f4-f26d-4d92-b221-cff884052dfa" alt="Breizorro Logo">
</p>

Breizorro is a flexible software program made to simplify image analysis tasks, including identifying emission islands and generating and modifying image masks, which are frequently used in radio interferometric imaging.


# Parameter definition

```
Usage: breizorro [OPTIONS]

Options:
  -r, --restored-image IMAGE      Restored image file from which to build the
                                  mask
  -m, --mask-image MASK           Input mask file(s). Either restored-image or
                                  mask-image must be specified.
  -t, --threshold float           Sigma threshold for masking (default = 6.5)
  -b, --boxsize int               Box size over which to compute stats
                                  (default = 50)
  --savenoise / --no-savenoise    Export noise image as FITS file
  --merge MASK(s)|REG(s)          Merge one or more masks or region files
  --subtract MASK(s)|REG(s)       Subtract one or more masks or region files
  -rc, --radial-cutoff RADIUS     Zero out mask beyond this radius in pixels
                                  from center (imitates beam attenuation)
  --number-islands / --no-number-islands
                                  Number the islands detected
  --remove-islands N|COORD        Remove islands from input mask (list by
                                  number or coordinates) e.g. --remove-islands
                                  1,18,20,20h10m13s:14d15m20s
  --ignore-missing-islands / --no-ignore-missing-islands
                                  Do not throw an error if an island specified
                                  by coordinates does not exist
  --extract-islands N|COORD       Extract islands from input mask (list by
                                  number or coordinates) e.g. --extract-
                                  islands 1,18,20,20h10m13s:14d15m20s
  --minimum-size int              Remove islands with areas fewer than or
                                  equal to the specified number of pixels
  --make-binary / --no-make-binary
                                  Replace all island numbers with 1
  --invert / --no-invert          Invert the mask
  --dilate R                      Apply dilation with a radius of R pixels
  --erode N                       Apply N iterations of erosion
  --fill-holes / --no-fill-holes  Fill holes (closed regions) in the mask
  --sum-peak float                Sum-to-peak ratio of flux islands to mask in
                                  original image
  -j, --ncpu int                  Number of processors to use for cataloguing
  --beam-size float               Average beam size in arcsec if missing in
                                  the image header
  -sf, --source-fitting str       Source fitting method for catalogue
                                  generation
  --gui / --no-gui                Open mask in bokeh html gui
  -o, --outfile File              The output mask image (default based on
                                  input name)
  --outcatalog File               Generate a catalogue based on the region
                                  mask
  --outregion File                Generate polygon regions from the mask
  --help                          Show this message and exit.
```


# Image Masking and Transformation

Breizorro uses the minimal filter to generate binary masks by replacing each pixel's value with the lowest value found in its near vicinity. This is incredibly effective at reducing noise and highlighting or smoothing particular regions of an image, such as the regions surrounding bright or small sources. Users can specify a window (or kernel) of a specific size that moves over the image as well as a sigma threshold for masking.

```
breizorro -t 6.5 -b 50 -r circinus-MFS-image.fits
```

![mypipelinerun_circinus_p3_3-MFS-image fits-mypipelinerun_circinus_p3_3-MFS-image mask fits-image-2024-09-11-10-35-43](https://github.com/user-attachments/assets/cfd2f918-340a-4148-96a2-c00ca41b33d0)


```
breizorro -r circinus-MFS-image.fits --sum-peak 500
```

![mypipelinerun_circinus_p3_3-MFS-image mask fits-mypipelinerun_circinus_p3_3-MFS-image mask fits-image-2024-09-11-13-23-14](https://github.com/user-attachments/assets/0ff50068-ec8a-42bf-8539-9f68f15a1ea9)

# Region Generation and Manipulation

Breizorro makes it easier to create and work with regions using image masks. It includes labelling, eliminating, extracting, and filtering regions (islands) based on user-specified criteria. Users can refine their regions of interest using techniques such as erosion, dilation, hole-filling, binary masking, and inversion.

```
breizorro -r circinus-MFS-image.fits --outregion circinus.reg
```

![mypipelinerun_circinus_p3_3-MFS-image fits-image-2024-09-11-10-38-15](https://github.com/user-attachments/assets/14f435e1-6234-4515-9597-c3002a644975)

```
breizorro -r circinus-MFS-image.fits --merge west.reg --dilate 1 --fill-holes
```

![mypipelinerun_circinus_p3_3-MFS-image mask fits-mypipelinerun_circinus_p3_3-MFS-image mask fits-image-2024-09-11-13-59-39](https://github.com/user-attachments/assets/2308c7b7-2ec0-4895-b93b-5d96f3d99337)

# Radial Masking

Breizorro includes radial cutoff capability for beam attenuation simulation.

## Radial Cutoff

Simulate beam attenuation by zeroing everything beyond a specified circular radius:

```
breizorro -r circinus-MFS-image.fits --radial-cutoff 200
```

The `--radial-cutoff` (or `-rc`) parameter accepts a radius in pixels from the centre of the image. This creates a circular mask that zeros out all pixels beyond the specified radius, which is useful for:

- Imitating the natural beam attenuation in radio observations
- Removing edge effects and artefacts
- Focusing analysis on the central, more sensitive region of observations

The radial cutoff is always applied **after** all other mask operations.

# Cataloguing and Visualisation

Breizorro enables catalogue generation from extracted regions, saving source properties to an ASCII/text file.
This is particularly useful for analysing fields dominated by point sources.
By efficiently parameterising and cataloguing compact sources, Breizorro enables rapid cross-matching.

```
breizorro -r circinus-MFS-image.fits --outcatalog deep2.txt
```

```
# processing fits image: circinus-MFS-image.fits
# mean beam size (arcsec): 37.97
# original image peak flux (Jy/beam): 0.7683289647102356
# noise out (µJy/beam): 736.96
# cutt-off flux  (mJy/beam): 4.79
# freq0 (Hz): 1419944335.9375
# number of sources detected: 35
#
#format: name ra_d dec_d i i_err i_peak i_peak_error emaj_s emin_s pa_d
src0 212.1114852668253 -65.52022230002818 0.01645 0.0017 0.01645 0.0017 0.0 0.0 0.0
src1 212.15838150612473 -65.44043978204705 0.00513 0.00055 0.00513 0.00055 0.0 0.0 0.0
src2 212.20192908332638 -65.39404261881141 0.00537 0.00059 0.00537 0.00059 0.0 0.0 0.0
src3 212.25942030914388 -65.08587770640244 0.00894 0.00112 0.00894 0.00112 0.0 0.0 0.0
src4 212.34142138447163 -65.28756781209889 0.06875 0.00699 0.06875 0.00699 0.0 0.0 0.0
src5 212.4017279869488 -65.28100679894399 0.00814 0.00102 0.00814 0.00102 0.0 0.0 0.0
src6 212.48217912947138 -65.72299674107823 0.05203 0.00531 0.05203 0.00531 0.0 0.0 0.0
src7 212.50820130212415 -64.88584945077669 0.00536 0.0006 0.00536 0.0006 0.0 0.0 0.0
src8 212.5702902342093 -65.65050370039532 0.00546 0.00059 0.00546 0.00059 0.0 0.0 0.0
src9 212.57060833616958 -65.44161961395122 0.0081 0.00085 0.0081 0.00085 0.0 0.0 0.0
src10 212.71633807796545 -65.42932839276426 0.00483 0.00049 0.00483 0.00049 0.0 0.0 0.0
src11 212.77548695536805 -65.2127221382807 0.00502 0.00053 0.00502 0.00053 0.0 0.0 0.0
src12 212.94697030800972 -65.46771419775158 0.0059 0.00063 0.0059 0.00063 0.0 0.0 0.0
src13 213.02987212660275 -65.82402586988974 0.0058 0.00066 0.0058 0.00066 0.0 0.0 0.0
src14 213.09132783845993 -65.83601937078244 0.00893 0.00094 0.00893 0.00094 0.0 0.0 0.0
src15 213.09865015952178 -64.84001483960485 0.00919 0.00099 0.00919 0.00099 0.0 0.0 0.0
src16 213.11251590492353 -65.02758772062043 0.00651 0.00085 0.00651 0.00085 0.0 0.0 0.0
src17 213.23985204575342 -65.03719936384542 0.00538 0.0006 0.00538 0.0006 0.0 0.0 0.0
src18 213.2413293070179 -65.01106379809042 0.00533 0.00063 0.00533 0.00063 0.0 0.0 0.0
src19 213.2899899919003 -65.34038197853626 1.56983 0.15704 0.76833 0.07696 238.62 74.32 -41.88
src20 213.30053772896596 -65.55087714062736 0.00831 0.001 0.00831 0.001 0.0 0.0 0.0
src21 213.30590547065324 -65.23402090775575 0.00511 0.00056 0.00511 0.00056 0.0 0.0 0.0
src22 213.3703871930737 -65.02801216160296 0.00534 0.00062 0.00534 0.00062 0.0 0.0 0.0
src23 213.40473320520746 -65.7032754189847 0.02391 0.00261 0.02255 0.00249 40.22 27.75 135.34
src24 213.48631194937263 -65.78499670808618 0.00506 0.00051 0.00506 0.00051 0.0 0.0 0.0
src25 213.5033368507256 -65.0479188751319 0.04865 0.00499 0.04865 0.00499 0.0 0.0 0.0
```


---

### 1. Source Fitting Methods Added

**New Parameter**: `--source-fitting` / `-sf`

Four methods are now available:
- **centroid** (DEFAULT) - Centre of mass, fastest, most robust
- **gaussian** - 2D Gaussian fitting, highest accuracy, can hang (auto-fallback)
- **moments** - Image moments, a fast and robust alternative
- **windowed** - Windowed centroid around peak, balanced

---

## IMPORTANT: Installation

### Catalogue features require additional dependencies!

```bash
# For catalogue generation (RECOMMENDED)
pip install breizorro[catalog]

# For ALL features (catalog + GUI)
pip install breizorro[all]

# With uv (faster)
uv pip install breizorro[catalog]
uv pip install breizorro[all]
```

### What's Included:

| Package | Includes | Use Case |
|---------|----------|----------|
| `breizorro` | Basic masking | Simple mask operations |
| `breizorro[catalog]` | + photutils, scikit-image | **Source cataloging (MOST USERS)** |
| `breizorro[gui]` | + bokeh | Interactive GUI |
| `breizorro[all]` | Everything above | Full feature set |

---

## Advanced Usage Examples

### Catalog Generation (Requires [catalog])

```bash
# Install first!
pip install breizorro[catalog]

# Centroid fitting (default, fastest, most robust)
breizorro -r image.fits --outcatalog sources.cat

# Gaussian fitting (high accuracy, may hang)
breizorro -r image.fits --outcatalog sources.cat --source-fitting gaussian

# Moments fitting (fast alternative)
breizorro -r image.fits --outcatalog sources.cat --source-fitting moments

# With regions output
breizorro -r image.fits --outcatalog sources.cat --outregion sources.rgn
```

---

## Source Fitting Method Comparison

| Method | Speed | Accuracy | Robustness | Hangs? | When to Use |
|--------|-------|----------|------------|--------|-------------|
| **centroid** | ⚡⚡⚡ | ✓✓ | ✓✓✓ | Never | **Default, general use** |
| gaussian | ⚡ | ✓✓✓ | ✓ | Sometimes | Point sources, high accuracy |
| moments | ⚡⚡⚡ | ✓✓ | ✓✓✓ | Never | Fast alternative |
| windowed | ⚡⚡ | ✓✓ | ✓✓ | Never | Bright compact sources |

---

# Contributors

Thank you to the people who have contributed to this project.

[![Contributors](https://contrib.rocks/image?repo=ratt-ru/breizorro)](https://github.com/ratt-ru/breizorro/graphs/contributors)
