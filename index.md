
---
layout: default
---

![Octocat](https://github.githubassets.com/images/icons/emoji/octocat.png)

Breizorro is a flexible software program made to simplify image analysis tasks, including identifying emission islands and generating and modifying picture masks, which are frequently used in radio interferometric imaging.

# Image Masking and Transformation

Breizorro uses the minimal filter to generate binary masks by replacing each pixel's value with the lowest value found in its near vicinity. This is incredibly effective at reducing noise and highlighting or smoothing particular regions of a picture, such as the regions surrounding bright or small sources. Users can specify a window (or kernel) of a specific size that moves over the image as well as a sigma threshold for masking.

# Region Generation and Manipulation

Breizorro makes it easier to create and work with regions using image masks. Labeling, eliminating, extracting, and filtering regions (islands) based on user-specified criteria are all included in this. Users can employ techniques including erosion, dilation, hole-filling, binary masking, and inversion to refine their regions of interest.

# Cataloguing and Visualization

Breizorro enables catalog generation from extracted regions, saving source properties to an ASCII/text
file.

```
breizorro -r circinus-MFS-image.fits --save-catalog circinus.txt
```
