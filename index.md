<p align="center">
  <img src="https://github.com/user-attachments/assets/478422f4-f26d-4d92-b221-cff884052dfa" alt="Breizorro Logo">
</p>

Breizorro is a flexible software program made to simplify image analysis tasks, including identifying emission islands and generating and modifying image masks, which are frequently used in radio interferometric imaging.


# Parameter definition

```
Usage: breizorro [OPTIONS]

Options:
  --restored-image PATH           Restored image file from which to build the
                                  mask
  --mask-image PATH               Input mask file(s). Either restored-image or
                                  mask-image must be specified.
  --threshold FLOAT               Sigma threshold for masking (default = 6.5)
  --boxsize INTEGER               Box size over which to compute stats
                                  (default = 50)
  --savenoise / --no-savenoise    Export noise image as FITS file
  --merge TEXT                    Merge one or more masks or region files
  --subtract TEXT                 Subtract one or more masks or region files
  --number-islands / --no-number-islands
                                  Number the islands detected
  --remove-islands TEXT           Remove islands from input mask (list by
                                  number or coordinates)
  --ignore-missing-islands / --no-ignore-missing-islands
                                  Do not throw an error if an island specified
                                  by coordinates does not exist
  --extract-islands TEXT          Extract islands from input mask (list by
                                  number or coordinates)
  --minimum-size INTEGER          Remove islands with areas fewer than or
                                  equal to the specified number of pixels
  --make-binary / --no-make-binary
                                  Replace all island numbers with 1
  --invert / --no-invert          Invert the mask
  --dilate INTEGER                Apply dilation with a radius of R pixels
  --erode INTEGER                 Apply N iterations of erosion
  --fill-holes / --no-fill-holes  Fill holes (closed regions) in the mask
  --sum-peak FLOAT                Sum-to-peak ratio of flux islands to mask in
                                  original image
  --ncpu INTEGER                  Number of processors to use for cataloging
  --beam-size FLOAT               Average beam size in arcsec if missing in
                                  the image header
  --gui / --no-gui                Open mask in bokeh html gui
  --outfile PATH                  Suffix for the mask image (default based on
                                  input name)
  --outcatalog PATH               Generate a catalog based on the region mask
  --outregion PATH                Generate polygon regions from the mask
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


# Cataloguing and Visualization

Breizorro enables catalogue generation from extracted regions, saving source properties to an ASCII/text file.
This is particularly useful for analyzing fields dominated by point sources.
By efficiently parameterizing and cataloguing compact sources, Breizorro enables rapid cross-matching.

```
breizorro -r deep2-MFS-image.fits --outcatalog deep2.txt
```

```
# processing fits image: deep2-MFS-image.fits
# mean beam size (arcsec): 6.31 
# original image peak flux (Jy/beam): 0.023107271641492844 
# noise out (µJy/beam): 51.62 
# cutt-off flux  (mJy/beam): 0.52 
# freq0 (Hz): 1049833007.8125 
# number of sources detected: 100 
#
#format: name ra_d dec_d i i_err emaj_s emin_s pa_d
src0 60.64975237181635 -79.86348735299585 0.00369 0.0001 11.66 0.0 120.96
src1 60.67140679629887 -80.36126947232378 7e-05 0.0001 0.0 0.0 0.0
src2 60.877334392876136 -79.76691180042988 0.00113 0.0001 8.25 0.0 104.04
src3 60.894387589282964 -79.73476060235502 0.00023 0.0001 0.0 0.0 0.0
src4 60.95196895741357 -79.68021884337942 0.00054 0.0001 0.0 0.0 0.0
src5 61.00657438220518 -80.17745581694626 0.00018 0.0001 0.0 0.0 0.0
src6 61.04136845645677 -80.48097816446368 0.00032 0.0001 0.0 0.0 0.0
src7 61.061502553348895 -79.93907167920088 0.01435 0.0001 21.54 0.0 111.8
src8 61.18650068905602 -80.40952881341589 0.00027 0.0001 0.0 0.0 0.0
src9 61.32594143681846 -79.91488530476678 7e-05 0.0001 0.0 0.0 0.0
src10 61.47462421089555 -79.96912396162651 9e-05 0.0001 0.0 0.0 0.0
src11 61.56419377531346 -80.18249455745902 0.0022 0.0001 10.77 0.0 158.2
```

# Contributors

Thank you to the people who have contributed to this project.

[![Contributors](https://contrib.rocks/image?repo=ratt-ru/breizorro)](https://github.com/ratt-ru/breizorro/graphs/contributors)
