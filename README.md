# Party

This folder holds the interesting implementation ideas
of particle physics entities.

## Four-vector

### Calculating pz with $eta$ ($theta$)

According to the documentation of `TLorentzVector` [$^1$](https://root.cern.ch/doc/master/TLorentzVector_8h_source.html#l00351)
$$p_{z} = p_{\text{T}}\sinh{\eta},$$
which might look strange at first shot.

This can be derived using the following relations
$$\eta = -\ln{\tan{\frac{\theta}{2}}},$$
$$p_{\text{T}} = p_{z}\tan{\theta}.$$

One can immediately get
$$p_{z} = \frac{p_{\text{T}}}{\tan{\big(2\arctan{(e^{-\eta})}\big)}}.$$
> The `Vector3` class still uses this formula [$^2$](https://root.cern.ch/doc/master/TVector3_8cxx_source.html#l00338)

Let $e^{-\eta} = \tan{\zeta}$, then
$$\frac{1}{\tan{\big(2\arctan{(e^{-\eta})}\big)}} = \frac{1}{\tan{2\zeta}}.$$

Remind trigonometric functions and substitute $eta$ back
$$\tan{2\zeta} = \frac{2\tan{\zeta}}{1-\tan^2{\zeta}} = \frac{2e^{-\eta}}{1-e^{-2\eta}} = \frac{2}{e^{\eta}-e^{-\eta}} = \frac{1}{\sinh{\eta}}.$$

Then the relation is proved $p_{z} = p_{\text{T}}\sinh{\eta}$.

### Dev Log

- reference: [Math::LorentzVector](https://root.cern.ch/doc/master/classROOT_1_1Math_1_1LorentzVector.html)
- delta_eta, delta_phi, delta_r(pseudorapidity: bool)
