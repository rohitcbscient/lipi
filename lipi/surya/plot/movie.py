from images2gif import writeGif
from PIL import Image
import os
import imageio

#file_names = sorted((fn for fn in os.listdir('.') if fn.endswith('.png')))[0:1000]
#['animationframa.png', 'animationframb.png', ...] "

def write(file_names,filename):
    images = [Image.open(fn) for fn in file_names]
    size = (2000,2000)
    for im in images:
        print im
        im.thumbnail(size, Image.ANTIALIAS)
    print writeGif.__doc__
    writeGif(filename, images, duration=0.001)


def write_imagero(file_names,output,frame):
    kwargs_write = {'fps':frame,'quantizer':'nq'}
    images = []
    for filename in file_names:
        print filename
        images.append(imageio.imread(filename))
    imageio.mimsave(output, images, 'GIF-FI', **kwargs_write)

