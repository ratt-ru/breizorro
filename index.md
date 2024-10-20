<p align="center">
  <img src="https://github.com/user-attachments/assets/478422f4-f26d-4d92-b221-cff884052dfa" alt="Breizorro Logo">
</p>

Breizorro is a flexible software program made to simplify image analysis tasks, including identifying emission islands and generating and modifying image masks, which are frequently used in radio interferometric imaging.


# Parameter definition

```
     breizorro [options] --restored-image restored_image

     optional arguments:
          -h, --help            show this help message and exit
          -r IMAGE, --restored-image IMAGE
                                Restored image file from which to build mask
          -m MASK, --mask-image MASK
                                Input mask file(s). Either --restored-image or --mask-
                                image must be specfied.
          -t THRESHOLD, --threshold THRESHOLD
                                Sigma threshold for masking (default = 6.5)
          -b BOXSIZE, --boxsize BOXSIZE
                                Box size over which to compute stats (default = 50)
          --savenoise           Enable to export noise image as FITS file (default=do
                                not save noise image)
          --merge MASK(s)|REG(s) [MASK(s)|REG(s) ...]
                                Merge in one or more masks or region files
          --subtract MASK(s)|REG(s) [MASK(s)|REG(s) ...]
                                Subract one or more masks or region files
          --number-islands      Number the islands detected (default=do not number
                                islands)
          --remove-islands N|COORD [N|COORD ...]
                                List of islands to remove from input mask. e.g.
                                --remove-islands 1 18 20 20h10m13s,14d15m20s
          --ignore-missing-islands
                                If an island specified by coordinates does not exist,
                                do not throw an error
          --extract-islands N|COORD [N|COORD ...]
                                List of islands to extract from input mask. e.g.
                                --extract-islands 1 18 20 20h10m13s,14d15m20s
          --minimum-size MINSIZE
                                Remove islands that have areas fewer than or equal to
                                the specified number of pixels
          --make-binary         Replace all island numbers with 1
          --invert              Invert the mask
          --dilate R            Apply dilation with a radius of R pixels
          --erode N             Apply N iterations of erosion
          --fill-holes          Fill holes (i.e. entirely closed regions) in mask
          --sum-peak SUM_PEAK   Sum to peak ratio of flux islands to mask in original
                                image.e.g. --sum-peak 100 will mask everything with a
                                ratio above 100
          -o OUTFILE, --outfile OUTFILE
                                Suffix for mask image (default based on input name
          --gui                 Open mask in gui.
```


# Image Masking and Transformation

Breizorro uses the minimal filter to generate binary masks by replacing each pixel's value with the lowest value found in its near vicinity. This is incredibly effective at reducing noise and highlighting or smoothing particular regions of a picture, such as the regions surrounding bright or small sources. Users can specify a window (or kernel) of a specific size that moves over the image as well as a sigma threshold for masking.

```
breizorro -t 6.5 -b 50 -r circinus-MFS-image.fits
```

![mypipelinerun_circinus_p3_3-MFS-image fits-mypipelinerun_circinus_p3_3-MFS-image mask fits-image-2024-09-11-10-35-43](https://github.com/user-attachments/assets/cfd2f918-340a-4148-96a2-c00ca41b33d0)


```
breizorro -r circinus-MFS-image.fits --sum-peak 500
```

![mypipelinerun_circinus_p3_3-MFS-image mask fits-mypipelinerun_circinus_p3_3-MFS-image mask fits-image-2024-09-11-13-23-14](https://github.com/user-attachments/assets/0ff50068-ec8a-42bf-8539-9f68f15a1ea9)

# Region Generation and Manipulation

Breizorro makes it easier to create and work with regions using image masks. It includes labeling, eliminating, extracting, and filtering regions (islands) based on user-specified criteria. Users can refine their regions of interest using techniques such as erosion, dilation, hole-filling, binary masking, and inversion.

```
breizorro -r circinus-MFS-image.fits --save-regions circinus.reg
```

![mypipelinerun_circinus_p3_3-MFS-image fits-image-2024-09-11-10-38-15](https://github.com/user-attachments/assets/14f435e1-6234-4515-9597-c3002a644975)

```
breizorro -r circinus-MFS-image.fits --merge west.reg --dilate 1 --fill-holes
```

![mypipelinerun_circinus_p3_3-MFS-image mask fits-mypipelinerun_circinus_p3_3-MFS-image mask fits-image-2024-09-11-13-59-39](https://github.com/user-attachments/assets/2308c7b7-2ec0-4895-b93b-5d96f3d99337)


# Cataloguing and Visualization

Breizorro enables catalog generation from extracted regions, saving source properties to an ASCII/text
file.

```
breizorro -r circinus-MFS-image.fits --save-catalog circinus.txt
```

```
# processing fits image: circinus-MFS-image.fits
# mean beam size (arcsec): 37.97
# original image peak flux (Jy/beam): 0.768315315246582
# noise out (μJy/beam): 737.33
# cutt-off flux (mJy/beam): 4.79
# freq0 (Hz): 1419944335.9375
# number of sources detected: 35
#
#format: name ra_d dec_d i i_err emaj_s emin_s pa_d
src0 -147.88904620760084 -65.52043870236054 0.00468 0.0001 50.04 0.0 177.71
src1 -147.84162114219356 -65.44040921224196 0.00047 0.0001 0.0 0.0 0.0
src2 -147.7980562104931 -65.39402001087845 0.00089 0.0001 0.0 0.0 0.0
src3 -147.74069959466678 -65.08592423258469 0.0082 0.0001 51.61 0.0 144.46
src4 -147.6583486798732 -65.28753961408073 0.09922 0.0001 86.58 0.0 -173.37
src5 -147.59829620974554 -65.28104296056243 0.00629 0.0001 47.07 0.0 167.74
src6 -147.51788583319714 -65.72316001301358 0.06639 0.0001 71.47 0.0 162.07
src7 -147.49179712006742 -64.88584736668952 0.00092 0.0001 0.0 0.0 0.0
src8 -147.42970139783463 -65.65050232604096 0.00069 0.0001 0.0 0.0 0.0
src9 -147.4293961973031 -65.44158144271519 0.00134 0.0001 0.0 0.0 0.0
src10 -147.28370646739054 -65.42936506202037 9e-05 0.0001 0.0 0.0 0.0
```

```
breizorro -r circinus-MFS-image.fits --save-catalog circinus.txt
```

![Screenshot 2024-09-11 100116](https://github.com/user-attachments/assets/79d3ece6-5d96-48a1-a06e-fbae2987d333)
