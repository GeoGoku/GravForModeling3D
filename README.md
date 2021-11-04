Abstract
Forward modeling is a requisite for many geophysical applications; for instance, forward modeling is used to compute the gravity anomalies in local- or regional-scale geological models. In this study, we develop a simple and efficient approach for the 3D interactive and fast modeling. The subsurface is divided into rectangular cells with constant density. Then, the cells are merged into a series of cell groups to reduce the storage overhead. Finally, the gravity anomalies at survey locations are computed using these cell groups. All the implementations presented in this paper are developed in MATLAB and executed in parallel. For synthetic models exceeding 107 cells with 4096 survey locations, the computation time of the proposed approach on a single GPU achieves an acceleration of 13 times compared to the “summation” approach on a single GPU, and 681 times compared to the “summation” approach on a single CPU core. The gravity anomalies calculated from the proposed approach and analytical solution are consistent, the RMS error is 1.4 × 10−12 mGal. For SEG/EAGE salt model, the proposed approach achieves a compression ratio of 855 for the model and takes 5.06 min to obtain gravity anomalies.
Keywords: Gravity; Forward modeling; Cell mergence; Parallel computing


Citation
Tao Chen, Guibin Zhang,
Forward modeling of gravity anomalies based on cell mergence and parallel computing,
Computers & Geosciences,
Volume 120,
2018,
Pages 1-9,
ISSN 0098-3004,
https://doi.org/10.1016/j.cageo.2018.07.007.


Please report any bug to geogoku@aliyun.com
