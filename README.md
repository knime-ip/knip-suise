Supervised Image Segmentation (experimental)
==========

The 'Supervised Image Segmentation'-plugin is a result of the Phd-Thesis from Martin Horn about the Active Segmentation of Images, to be published in 2014 at the University of Konstanz.

The thesis proposes a general Active Segmentation Framework that encompasses three main steps, the pixel model, segmentation, and the segment model:

1. Pixel Model: All pixels are classified beeing part of a certain region or not; or lying on the objects contour (boundary) or not.

2. Segmentation: The pixel classification gives a canonical representation of the orginal image that now allows the application of standard segmentation algorithms (e.g. watershed and region merging, connected component analysis, etc.). Instead of a (planar) segmentation this steps returns a suffienctly large set of segments (possibly overlapping) which is beeing processed in the next step.

3. Segment Model: The segments derived in the previous step are finally classified and weighted to the effect that desired segments receive a high probability of beeing positive. This allows one to compose the final segmentation by only considering positive segments of high probability and neglecting segments that overlap with other segments of higher probability.

The plugin provides a set of nodes that allows the realization of the aforementioned steps (e.g. the pixel-wise classification, either of boundaries or regions).
This [workflow] (http://tech.knime.org/files/knimeip/workflows/Example_Supervised Image Segmentation.zip) exemplary demonstrates the usage.
