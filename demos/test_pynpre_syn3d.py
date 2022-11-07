import pynpre as npre
import numpy as np
din=np.random.randn(1001)
dout=npre.npre3d(din)

err=np.linalg.norm(din-dout)
print(err)


