 from astropy.io import fits
>>> imax = fits.open('imax_find_sun_input.fits')[0].data
>>> imax.shape
(4, 5, 936, 936)
>>> 
>>> imax_i = imax[0,4,:,:]
>>> imax_i.shape
(936, 936)
>>> from matplotlib import pyplot as plt
>>> plt.imshow(imax_i, cmap='gray')
<matplotlib.image.AxesImage object at 0x7fceca9c7b70>
>>> plt.show()
>>> from scipy.fftpack import fftn
>>> i_fourier = fftn(imax_i)
Traceback (most recent call last):
  File "/home/prabhu/anaconda3/lib/python3.6/site-packages/scipy/fftpack/basic.py", line 630, in _raw_fftn_dispatch
    work_function = _DTYPE_TO_FFTN[tmp.dtype]
KeyError: dtype('>f4')

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/home/prabhu/anaconda3/lib/python3.6/site-packages/scipy/fftpack/basic.py", line 623, in fftn
    return _raw_fftn_dispatch(x, shape, axes, overwrite_x, 1)
  File "/home/prabhu/anaconda3/lib/python3.6/site-packages/scipy/fftpack/basic.py", line 632, in _raw_fftn_dispatch
    raise ValueError("type %s is not supported" % tmp.dtype)
ValueError: type >f4 is not supported
>>> imax_i.shape
(936, 936)
>>> imax_i.type
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
AttributeError: 'numpy.ndarray' object has no attribute 'type'
>>> imax_i.mean()
3484.7219
>>> 
>>> 
>>> from scipy import fftpack as fft
>>> fft.fftn(imax_i)
Traceback (most recent call last):
  File "/home/prabhu/anaconda3/lib/python3.6/site-packages/scipy/fftpack/basic.py", line 630, in _raw_fftn_dispatch
    work_function = _DTYPE_TO_FFTN[tmp.dtype]
KeyError: dtype('>f4')

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/home/prabhu/anaconda3/lib/python3.6/site-packages/scipy/fftpack/basic.py", line 623, in fftn
    return _raw_fftn_dispatch(x, shape, axes, overwrite_x, 1)
  File "/home/prabhu/anaconda3/lib/python3.6/site-packages/scipy/fftpack/basic.py", line 632, in _raw_fftn_dispatch
    raise ValueError("type %s is not supported" % tmp.dtype)
ValueError: type >f4 is not supported
>>> fft.fft2(imax_i)
Traceback (most recent call last):
  File "/home/prabhu/anaconda3/lib/python3.6/site-packages/scipy/fftpack/basic.py", line 630, in _raw_fftn_dispatch
    work_function = _DTYPE_TO_FFTN[tmp.dtype]
KeyError: dtype('>f4')

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/home/prabhu/anaconda3/lib/python3.6/site-packages/scipy/fftpack/basic.py", line 675, in fft2
    return fftn(x,shape,axes,overwrite_x)
  File "/home/prabhu/anaconda3/lib/python3.6/site-packages/scipy/fftpack/basic.py", line 623, in fftn
    return _raw_fftn_dispatch(x, shape, axes, overwrite_x, 1)
  File "/home/prabhu/anaconda3/lib/python3.6/site-packages/scipy/fftpack/basic.py", line 632, in _raw_fftn_dispatch
    raise ValueError("type %s is not supported" % tmp.dtype)
ValueError: type >f4 is not supported
>>> from scipy.fftpack import fftn, ifftn
>>> 
>>> 
>>> 
>>> arr = np.zeros((100,100))
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
NameError: name 'np' is not defined
>>> import numpy as np
>>> arr = np.zeros((100,100))
>>> 
>>> 
>>> arr.shape
(100, 100)
>>> 
>>> 
>>> fftn(arr)
array([[ 0.+0.j,  0.+0.j,  0.+0.j, ...,  0.+0.j,  0.+0.j,  0.+0.j],
       [ 0.+0.j,  0.+0.j,  0.+0.j, ...,  0.+0.j,  0.+0.j,  0.+0.j],
       [ 0.+0.j,  0.+0.j,  0.+0.j, ...,  0.+0.j,  0.+0.j,  0.+0.j],
       ..., 
       [ 0.+0.j,  0.+0.j,  0.+0.j, ...,  0.+0.j,  0.+0.j,  0.+0.j],
       [ 0.+0.j,  0.+0.j,  0.+0.j, ...,  0.+0.j,  0.+0.j,  0.+0.j],
       [ 0.+0.j,  0.+0.j,  0.+0.j, ...,  0.+0.j,  0.+0.j,  0.+0.j]])
>>> imax_i.shape
(936, 936)
>>> dir(imax_i)
['T', '__abs__', '__add__', '__and__', '__array__', '__array_finalize__', '__array_interface__', '__array_prepare__', '__array_priority__', '__array_struct__', '__array_wrap__', '__bool__', '__class__', '__complex__', '__contains__', '__copy__', '__deepcopy__', '__delattr__', '__delitem__', '__dir__', '__divmod__', '__doc__', '__eq__', '__float__', '__floordiv__', '__format__', '__ge__', '__getattribute__', '__getitem__', '__gt__', '__hash__', '__iadd__', '__iand__', '__ifloordiv__', '__ilshift__', '__imatmul__', '__imod__', '__imul__', '__index__', '__init__', '__init_subclass__', '__int__', '__invert__', '__ior__', '__ipow__', '__irshift__', '__isub__', '__iter__', '__itruediv__', '__ixor__', '__le__', '__len__', '__lshift__', '__lt__', '__matmul__', '__mod__', '__mul__', '__ne__', '__neg__', '__new__', '__or__', '__pos__', '__pow__', '__radd__', '__rand__', '__rdivmod__', '__reduce__', '__reduce_ex__', '__repr__', '__rfloordiv__', '__rlshift__', '__rmatmul__', '__rmod__', '__rmul__', '__ror__', '__rpow__', '__rrshift__', '__rshift__', '__rsub__', '__rtruediv__', '__rxor__', '__setattr__', '__setitem__', '__setstate__', '__sizeof__', '__str__', '__sub__', '__subclasshook__', '__truediv__', '__xor__', 'all', 'any', 'argmax', 'argmin', 'argpartition', 'argsort', 'astype', 'base', 'byteswap', 'choose', 'clip', 'compress', 'conj', 'conjugate', 'copy', 'ctypes', 'cumprod', 'cumsum', 'data', 'diagonal', 'dot', 'dtype', 'dump', 'dumps', 'fill', 'flags', 'flat', 'flatten', 'getfield', 'imag', 'item', 'itemset', 'itemsize', 'max', 'mean', 'min', 'nbytes', 'ndim', 'newbyteorder', 'nonzero', 'partition', 'prod', 'ptp', 'put', 'ravel', 'real', 'repeat', 'reshape', 'resize', 'round', 'searchsorted', 'setfield', 'setflags', 'shape', 'size', 'sort', 'squeeze', 'std', 'strides', 'sum', 'swapaxes', 'take', 'tobytes', 'tofile', 'tolist', 'tostring', 'trace', 'transpose', 'var', 'view']
>>> dir(arr)
['T', '__abs__', '__add__', '__and__', '__array__', '__array_finalize__', '__array_interface__', '__array_prepare__', '__array_priority__', '__array_struct__', '__array_wrap__', '__bool__', '__class__', '__complex__', '__contains__', '__copy__', '__deepcopy__', '__delattr__', '__delitem__', '__dir__', '__divmod__', '__doc__', '__eq__', '__float__', '__floordiv__', '__format__', '__ge__', '__getattribute__', '__getitem__', '__gt__', '__hash__', '__iadd__', '__iand__', '__ifloordiv__', '__ilshift__', '__imatmul__', '__imod__', '__imul__', '__index__', '__init__', '__init_subclass__', '__int__', '__invert__', '__ior__', '__ipow__', '__irshift__', '__isub__', '__iter__', '__itruediv__', '__ixor__', '__le__', '__len__', '__lshift__', '__lt__', '__matmul__', '__mod__', '__mul__', '__ne__', '__neg__', '__new__', '__or__', '__pos__', '__pow__', '__radd__', '__rand__', '__rdivmod__', '__reduce__', '__reduce_ex__', '__repr__', '__rfloordiv__', '__rlshift__', '__rmatmul__', '__rmod__', '__rmul__', '__ror__', '__rpow__', '__rrshift__', '__rshift__', '__rsub__', '__rtruediv__', '__rxor__', '__setattr__', '__setitem__', '__setstate__', '__sizeof__', '__str__', '__sub__', '__subclasshook__', '__truediv__', '__xor__', 'all', 'any', 'argmax', 'argmin', 'argpartition', 'argsort', 'astype', 'base', 'byteswap', 'choose', 'clip', 'compress', 'conj', 'conjugate', 'copy', 'ctypes', 'cumprod', 'cumsum', 'data', 'diagonal', 'dot', 'dtype', 'dump', 'dumps', 'fill', 'flags', 'flat', 'flatten', 'getfield', 'imag', 'item', 'itemset', 'itemsize', 'max', 'mean', 'min', 'nbytes', 'ndim', 'newbyteorder', 'nonzero', 'partition', 'prod', 'ptp', 'put', 'ravel', 'real', 'repeat', 'reshape', 'resize', 'round', 'searchsorted', 'setfield', 'setflags', 'shape', 'size', 'sort', 'squeeze', 'std', 'strides', 'sum', 'swapaxes', 'take', 'tobytes', 'tofile', 'tolist', 'tostring', 'trace', 'transpose', 'var', 'view']
>>> imax_i = np.asarray(imax_i)
>>> fftn(imax_i)
Traceback (most recent call last):
  File "/home/prabhu/anaconda3/lib/python3.6/site-packages/scipy/fftpack/basic.py", line 630, in _raw_fftn_dispatch
    work_function = _DTYPE_TO_FFTN[tmp.dtype]
KeyError: dtype('>f4')

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/home/prabhu/anaconda3/lib/python3.6/site-packages/scipy/fftpack/basic.py", line 623, in fftn
    return _raw_fftn_dispatch(x, shape, axes, overwrite_x, 1)
  File "/home/prabhu/anaconda3/lib/python3.6/site-packages/scipy/fftpack/basic.py", line 632, in _raw_fftn_dispatch
    raise ValueError("type %s is not supported" % tmp.dtype)
ValueError: type >f4 is not supported
>>> imax_i = imax_i.astype('double')
>>> fftn(imax_i)
array([[  3.05295120e+09 +0.00000000e+00j,
         -4.45528321e+06 -3.97743810e+07j,
         -2.30864281e+06 -2.24633673e+07j, ...,
          8.97037882e+05 +7.55868541e+06j,
         -2.30864281e+06 +2.24633673e+07j,
         -4.45528321e+06 +3.97743810e+07j],
       [ -2.74862921e+08 +4.59451277e+08j,
          4.54499692e+05 -4.50381246e+07j,
          1.51455374e+07 -4.92938199e+07j, ...,
          3.17040833e+07 -3.94271849e+07j,
          4.03006463e+07 -4.20544416e+07j,
          5.80681573e+07 -4.04828348e+07j],
       [ -3.45686068e+07 +3.41872591e+08j,
         -8.42641488e+06 -1.25036175e+07j,
          1.38791061e+06 -2.17912880e+07j, ...,
          8.15512838e+06 -3.58171974e+07j,
          5.38385485e+06 -4.43471742e+07j,
          2.13552506e+07 -5.87010622e+07j],
       ..., 
       [  9.24529851e+07 -1.94898063e+08j,
         -2.13681331e+07 +4.82335392e+07j,
         -1.86307292e+07 +3.30545584e+07j, ...,
         -2.89143665e+06 +9.52290443e+06j,
         -4.23619958e+06 +1.01733692e+07j,
          6.22079190e+05 +3.52366233e+06j],
       [ -3.45686068e+07 -3.41872591e+08j,
          2.13552506e+07 +5.87010622e+07j,
          5.38385485e+06 +4.43471742e+07j, ...,
          2.73515014e+05 +2.50096257e+07j,
          1.38791061e+06 +2.17912880e+07j,
         -8.42641488e+06 +1.25036175e+07j],
       [ -2.74862921e+08 -4.59451277e+08j,
          5.80681573e+07 +4.04828348e+07j,
          4.03006463e+07 +4.20544416e+07j, ...,
          1.83588547e+07 +4.34432978e+07j,
          1.51455374e+07 +4.92938199e+07j,
          4.54499692e+05 +4.50381246e+07j]])
>>> i_fourier = fftn(imax_i)
>>> i_fourier.shape
(936, 936)
>>> plt.imshow(i_fourier.real, cmap='gray')
<matplotlib.image.AxesImage object at 0x7fceca5b66a0>
>>> plt.show()
>>> plt.imshow(i_fourier.real, cmap='gray', vmin=-10, vmax=10)
<matplotlib.image.AxesImage object at 0x7fced668d630>
>>> plt.show()
>>> plt.imshow(i_fourier.real, cmap='gray', vmin=-30, vmax=30)
<matplotlib.image.AxesImage object at 0x7fceca7d69b0>
>>> plt.show()
>>> plt.imshow(i_fourier.real, cmap='gray', vmin=-100, vmax=100)
<matplotlib.image.AxesImage object at 0x7fceca47b3c8>
>>> plt.show()
>>> plt.imshow(i_fourier, cmap='gray', vmin=-100, vmax=100)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/home/prabhu/anaconda3/lib/python3.6/site-packages/matplotlib/pyplot.py", line 3157, in imshow
    **kwargs)
  File "/home/prabhu/anaconda3/lib/python3.6/site-packages/matplotlib/__init__.py", line 1898, in inner
    return func(ax, *args, **kwargs)
  File "/home/prabhu/anaconda3/lib/python3.6/site-packages/matplotlib/axes/_axes.py", line 5124, in imshow
    im.set_data(X)
  File "/home/prabhu/anaconda3/lib/python3.6/site-packages/matplotlib/image.py", line 596, in set_data
    raise TypeError("Image data can not convert to float")
TypeError: Image data can not convert to float
>>> plt.imshow(i_fourier.imag, cmap='gray', vmin=-100, vmax=100)
<matplotlib.image.AxesImage object at 0x7fceca4008d0>
>>> plt.show()
>>> 
>>> 
>>> plt.imshow(i_fourier.real, cmap='gray', vmin=-100, vmax=100)
<matplotlib.image.AxesImage object at 0x7fceca395c18>
>>> plt.clf()
>>> plt.imshow(i_fourier.real, cmap='gray', vmin=-200, vmax=200)
<matplotlib.image.AxesImage object at 0x7fcec8b17940>
>>> plt.show()
>>> plt.imshow(imax_i, cmap='gray')
<matplotlib.image.AxesImage object at 0x7fcec8459e80>
>>> plt.show()
>>> i_fourier = imax_i[400:800,200:600]
>>> imax_crop = imax_i[400:800,200:600]
>>> plt.imshow(imax_crop, cmap='gray')
<matplotlib.image.AxesImage object at 0x7fcec8364240>
>>> plt.show()
>>> fou_crop = fftn(imax_crop)
>>> plt.imshow(fou_crop.real, cmap='gray', vmin=-100,vmax=100)
<matplotlib.image.AxesImage object at 0x7fcec81ba080>
>>> plt.show()
>>> plt.imshow(fou_crop.real, cmap='gray')
<matplotlib.image.AxesImage object at 0x7fcec80d05f8>
>>> plt.show()
>>> plt.imshow(fou_crop.real, cmap='gray', vmin=-500,vmax=500)
<matplotlib.image.AxesImage object at 0x7fcec8077240>
>>> fou_crop.mean()
(3861.56689453125-4.6566128730773928e-14j)
>>> np.std(fou_crop.real)
1708928.9097907722
>>> np.mean(fou_crop.real)
3861.5668945312491
>>> plt.imshow(fou_crop.real, cmap='gray', vmin=-50000,vmax=50000)
<matplotlib.image.AxesImage object at 0x7fcec8077828>
>>> plt.show()
>>> power = np.abs(fou_crop)
>>> fou_crop.shape
(400, 400)
>>> plt.imshow(power)
<matplotlib.image.AxesImage object at 0x7fcec812ce10>
>>> power.mean()
31121.731802815866
>>> power.std()
1712433.617166939
>>> plt.clf()
>>> plt.imshow(power, cmap='gray', vmin=-100000,vmax=100000)
<matplotlib.image.AxesImage object at 0x7fcec818d7f0>
>>> plt.show()
>>> ymean = fou_crop[:,150:250]
>>> plt.imshow(ymean, cmap='gray', vmin=-100000, vmax=100000)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/home/prabhu/anaconda3/lib/python3.6/site-packages/matplotlib/pyplot.py", line 3157, in imshow
    **kwargs)
  File "/home/prabhu/anaconda3/lib/python3.6/site-packages/matplotlib/__init__.py", line 1898, in inner
    return func(ax, *args, **kwargs)
  File "/home/prabhu/anaconda3/lib/python3.6/site-packages/matplotlib/axes/_axes.py", line 5124, in imshow
    im.set_data(X)
  File "/home/prabhu/anaconda3/lib/python3.6/site-packages/matplotlib/image.py", line 596, in set_data
    raise TypeError("Image data can not convert to float")
TypeError: Image data can not convert to float
>>> plt.imshow(ymean.real, cmap='gray', vmin=-100000, vmax=100000)
<matplotlib.image.AxesImage object at 0x7fceca400860>
>>> plt.show()
>>> ymean.shape
(400, 100)
>>> ym = np.mean(ymean, axis=0)
>>> ym.shape
(100,)
>>> ym = np.mean(ymean, axis=1)
>>> ym.shape
(400,)
>>> 
>>> 
>>> plt.plot(ym)
/home/prabhu/anaconda3/lib/python3.6/site-packages/numpy/core/numeric.py:531: ComplexWarning: Casting complex values to real discards the imaginary part
  return array(a, dtype, copy=False, order=order)
[<matplotlib.lines.Line2D object at 0x7fceca5ae978>]
>>> plt.plot(ym.real)
[<matplotlib.lines.Line2D object at 0x7fceca5ae3c8>]
>>> ym.shape
(400,)
>>> plt.show()
>>> fou_crop.real /= ym.real
>>> xmean = fou_crop[150:250,:]
>>> xm = np.mean(xmean, axis=1)
>>> xm.shape
(100,)
>>> xm = np.mean(xmean, axis=0)
>>> xm.shape
(400,)
>>> 
>>> 
>>> fou_crop.real /= xm.real
>>> plt.imshow(fou_crop.real, cmap='gray', vmin=-100000, vmax=100000)
<matplotlib.image.AxesImage object at 0x7fced668d198>
>>> plt.show()
>>> fou_crop = fftn(imax_crop)
>>> plt.imshow(fou_crop.real, cmap='gray', vmin=-100000, vmax=100000)
<matplotlib.image.AxesImage object at 0x7fcec2cccc50>
>>> plt.show9)
  File "<stdin>", line 1
    plt.show9)
             ^
SyntaxError: invalid syntax
>>> plt.show()
>>> x = fou_crop[200,:]
>>> plt.plot(x.real)
[<matplotlib.lines.Line2D object at 0x7fcec8190780>]
>>> plt.show()
>>> fou2 = fou_crop.real/x.real
>>> fou2.shape
(400, 400)
>>> 
>>> 
>>> plt.imshow(fou2, cmap='gray', vmin=-100000, vmax=100000)
<matplotlib.image.AxesImage object at 0x7fcec8af18d0>
>>> plt.show()
>>> plt.imshow(fou2, cmap='gray', vmin=-100, vmax=100)
<matplotlib.image.AxesImage object at 0x7fcec837b668>
>>> plt.show()
>>> xmean = fou_crop[150:250,:]
>>> xmean.real = np.mean(xmean.real, axis=0)
>>> xmean.real.shape
(100, 400)
>>> xmean.shape
(100, 400)
>>> 
>>> 
>>> xm = xmean.real
>>> xm.shape
(100, 400)
>>> 
>>> 
>>> xm = np.mean(xm, axis=0)
>>> xm.shape
(400,)
>>> plt.plot(xm)
[<matplotlib.lines.Line2D object at 0x7fcec845ef60>]
>>> plt.show()
>>> fou2 = fou_crop.real/xm
>>> plt.imshow(fou2, cmap='gray', vmin=-100, vmax=100)
<matplotlib.image.AxesImage object at 0x7fcec36ed0b8>
>>> plt.show()
>>> xmean = fou_crop[:,:].real
>>> xmean.shape
(400, 400)
>>> xm = np.mean(xmean, axis=0)
>>> xm.shape
(400,)
>>> fou2 = fou_crop.real/xm
>>> plt.imshow(fou2, cmap='gray', vmin=-100, vmax=100)
<matplotlib.image.AxesImage object at 0x7fceca359668>
>>> plt.show()
>>> plt.show()
>>> plt.imshow(fou2, cmap='gray', vmin=-100, vmax=100)
<matplotlib.image.AxesImage object at 0x7fceca59c588>
>>> plt.show()
>>> from scipy.fftpack import ifftn
>>> image = ifftn(fou2)
>>> plt.imshow(image,cmap='gray')
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/home/prabhu/anaconda3/lib/python3.6/site-packages/matplotlib/pyplot.py", line 3157, in imshow
    **kwargs)
  File "/home/prabhu/anaconda3/lib/python3.6/site-packages/matplotlib/__init__.py", line 1898, in inner
    return func(ax, *args, **kwargs)
  File "/home/prabhu/anaconda3/lib/python3.6/site-packages/matplotlib/axes/_axes.py", line 5124, in imshow
    im.set_data(X)
  File "/home/prabhu/anaconda3/lib/python3.6/site-packages/matplotlib/image.py", line 596, in set_data
    raise TypeError("Image data can not convert to float")
TypeError: Image data can not convert to float
>>> type(image)
<class 'numpy.ndarray'>
>>> image.shape
(400, 400)
>>> plt.imshow(image.real,cmap='gray')
<matplotlib.image.AxesImage object at 0x7fceca315f60>
>>> plt.show()
