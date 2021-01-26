## Multi-Modal-Image-Metric
The code of the manuscript "A Multi-Modal Edge Consistency Metric based on Regression Robustness of Truncated SVD".


## Function
`Func_ED.m` in the file "Funcs" is our metric.

`metricComparison_spectral.m` produces the experimental results of metric comparison on multi-spectral images.

`metricComparison_stereo.m` produces the experimental results of metric comparison on intensity/desparity images.

`metricComparison_spectral_buildSample.m` builds the test samples based on Multispectral Image Database of Columbia University [1].

`metricComparison_stereo_buildSample.m` builds the test samples based on Middlebury Stereo Dataset 2014 [2].

`stereoExample_EC.m` gives an example of using our metric for multi-modal stereo.

`stereoExample_DASIY.m` gives an example of using DASIY [3] for multi-modal stereo.

`stereoExample_DASC.m` gives an example of using DASC [4] for multi-modal stereo.

`stereoExample_L2net.m` gives an example of using L2-Net [5] for multi-modal stereo.

`build_DASC.m` is used to pre-calculate the DASC dense descriptors, which is too large to upload. 
We leave empty cookies in "Fig_metricComparison", thus in default, 
the metricComparison scripts will ignore DASC.
If you want to reproduce our results of DASC, you need to install DASC toolbox (https://github.com/seungryong/DASC) and run this script to build
"DASC spectral data.mat", "DASC stereo data.mat", and "DASC stereo example.mat".
You can find detailed instructions in this file.

`build_L2.m` is used to pre-calculate the L2-Net dense descriptors, which is too large to upload. 
We leave empty cookies in "Fig_metricComparison", thus in default, 
the metricComparison scripts will ignore L2-Net.
If you want to reproduce our results of L2-Net, you need to install matconvnet and  L2-Net toolbox (https://github.com/yuruntian/L2-Net), then run this script to build
"L2Net spectral data.mat", "L2Net stereo data.mat", and "L2Net stereo example.mat".
It will take a very long time (about 7 hours) if GPU is unavailable.
You can find detailed instructions in this file.

## References
[1] F. Yasuma, T. Mitsunaga, D. Iso, and S. Nayar, “Generalized assorted pixel camera: Post-capture control of resolution, dynamic range and spectrum,” Tech. Rep., Nov 2008.

[2] D. Scharstein, H. Hirschmu ̈ller, Y. Kitajima, G. Krathwohl, N. Nesˇic ́, X. Wang, and P. Westling, “High-resolution stereo datasets with subpixel-accurate ground truth,” in Proc. German Conference on Pattern Recognition, 2014, pp. 31–42.

[3] E. Tola, V. Lepetit, and P. Fua, “Daisy: An efficient dense descriptor applied to wide-baseline stereo,” IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 32, no. 5, pp. 815–830, 2009.

[4] S.Kim,D.Min,B.Ham,S.Ryu,M.N.Do,andK.Sohn,“Dasc:Dense adaptive self-correlation descriptor for multi-modal and multi-spectral correspondence,” in Proc. IEEE Conference on Computer Vision and Pattern Recognition, 2015, pp. 2103–2112.

[5] Y. Tian, B. Fan, and F. Wu, “L2-net: Deep learning of discriminative patch descriptor in euclidean space,” in Proc. IEEE Conference on Computer Vision and Pattern Recognition, 2017, pp. 661–669.
