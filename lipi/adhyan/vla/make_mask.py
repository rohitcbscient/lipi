

filename='image.im'
msname='20130423T2030-2050.L.50ms.selfcal.ms'
im.open('fantasticvis.ms')
im.drawmask(image=filename,mask='blank3.mask')
im.close()
