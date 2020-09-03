# Deep Generative Model of Phase Transformation in Layered Materials

[P. Rajak, *et al.*, *Phys. Rev. B 100*, 014108: 1-7 (2019)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.100.014108)

<img src="Highlight_Figure.png" width="545.6" height="192.8" align="right">

Optical and electrical properties of two-dimensional (2D) layered materials can be tuned by mechanical straining, which induces transformations from semiconducting to metallic phases. Here, we have developed a generative model of the phase transformation pathways in MoWSe<sub>2</sub> hetero-structure during  dyanamic fracture using variational autoencoder (VAE). The training dataset for VAE is generated using molecular dynamics (MD) simulation for different system sizes and strain rate. After training, VAE correctly identifies transformation pathways connecting the semiconducting (2H) and metallic (1T) phases via novel intermediate structures called _&alpha;_ and _&beta;_, which is also observed experimently.  Further, the stability of the structures synthesized from VAE are validated by quantum simulations based on density functional theory.</br>
This repository includes (1) C code to construct training dataset for VAE and Conditional VAE (CVAE) in numpy format from the MD simulaition trajectories of MoWSe<sub>2</sub>  fracture, (2) A sample fracture frame form MD simulation and (3) Ipython notebook for the training and anlysis of the VAE and CVAE models.
