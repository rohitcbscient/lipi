import numpy as np
from cora.signal import corr21cm
from cora.foreground import galaxy
import healpy
import matplotlib.pyplot as plt

cr = galaxy.FullSkySynchrotron()
aps1 = cr.angular_powerspectrum(np.arange(1000), 800.0, 800.0)
assert len(aps1) == 1000
assert np.allclose(aps1.sum(), 75.47681191093129, rtol=1e-7)  # Calculated for commit 02f4d1cd3f402d
fa = np.linspace(400.0, 800.0, 64)
aps2 = cr.angular_powerspectrum(np.arange(1000)[:, None, None], fa[None, :, None], fa[None, None, :])
assert aps2.shape == (1000, 64, 64)

gpolsky=cr.getpolsky()
gsky=cr.getsky()

healpy.visufunc.mollview(gsky[0])
plt.show()



