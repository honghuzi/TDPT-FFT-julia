# TDPT-xuv+ir
TDPT is simple

<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;P(E_1,&space;E_2)&space;=&space;1/2(\frac{c}{4\pi^2})^2|\sqrt{\frac{\sigma^{He}(E_1)\sigma^{He^&plus;}(E_2)}{\omega_{ai}\omega_{fa}}}K(E_a)&space;&plus;&space;\sqrt{\frac{\sigma^{He}(E_{2})\sigma^{He^&plus;}(E_1)}{\omega_{ai}\omega_{fa}}}K(E_b)|^2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;P(E_1,&space;E_2)&space;=&space;1/2(\frac{c}{4\pi^2})^2|\sqrt{\frac{\sigma^{He}(E_1)\sigma^{He^&plus;}(E_2)}{\omega_{ai}\omega_{fa}}}K(E_a)&space;&plus;&space;\sqrt{\frac{\sigma^{He}(E_{2})\sigma^{He^&plus;}(E_1)}{\omega_{ai}\omega_{fa}}}K(E_b)|^2" title="P(E_1, E_2) = 1/2(\frac{c}{4\pi^2})^2|\sqrt{\frac{\sigma^{He}(E_1)\sigma^{He^+}(E_2)}{\omega_{ai}\omega_{fa}}}K(E_a) + \sqrt{\frac{\sigma^{He}(E_{2})\sigma^{He^+}(E_1)}{\omega_{ai}\omega_{fa}}}K(E_b)|^2" /></a>

Here 

<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;K(E_a)=\int_{-\infty}^{\infty}\int_{-\infty}^{\tau_1}F(\tau_1)F(\tau_2)e^{i((Ea-Ei)*\tau_1&plus;(Ef-Ea)*\tau_2)}d\tau_1d\tau_2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;K(E_a)=\int_{-\infty}^{\infty}\int_{-\infty}^{\tau_1}F(\tau_1)F(\tau_2)e^{i((Ea-Ei)*\tau_1&plus;(Ef-Ea)*\tau_2)}d\tau_1d\tau_2" title="K(E_a)=\int_{-\infty}^{\infty}\int_{-\infty}^{\tau_1}F(\tau_1)F(\tau_2)e^{i((Ea-Ei)*\tau_1+(Ef-Ea)*\tau_2)}d\tau_1d\tau_2" /></a>

This double integral can be convert to a FFT via

<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;G(\tau_1,&space;\tau_2)&space;=&space;F(\tau_1)F(\tau_2)\times(\tau_1>\tau_2)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;G(\tau_1,&space;\tau_2)&space;=&space;F(\tau_1)F(\tau_2)\times(\tau_1>\tau_2)" title="G(\tau_1, \tau_2) = F(\tau_1)F(\tau_2)\times(\tau_1>\tau_2)" /></a>


The ac-stark shift is included both for xuv and ir.
