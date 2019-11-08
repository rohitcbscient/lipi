from images2gif import writeGif
from PIL import Image
import os


file_names = sorted((fn for fn in os.listdir('.') if fn.endswith('.png')))[0:1000]
#['animationframa.png', 'animationframb.png', ...] "

images = [Image.open(fn) for fn in file_names]

size = (2000,2000)
for im in images:
    im.thumbnail(size, Image.ANTIALIAS)

print writeGif.__doc__

filename = "20180421_240MHz.GIF"
writeGif(filename, images, duration=0.001)
