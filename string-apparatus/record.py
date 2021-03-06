#!/usr/bin/env python
import os
import datetime
now = datetime.datetime.now()
ext = 'mpg'

#size = 'hd1080'
#size = 'hd720'
size = 'xga'

outfile = now.strftime('screencaps/screencap-%d%h%YT%H%M%S.' + ext)
#cmd = 'ffmpeg -f x11grab -s ' + size + ' -r 30 -i :0.0  ' + outfile
#cmd = 'ffmpeg -f x11grab -s ' + size + ' -r 30 -i :0.0  -b 10000000 ' + outfile
#cmd = 'ffmpeg -f x11grab -s ' + size + ' -r 30 -i :0.0  -target ntsc-dvd ' + outfile

#audio_options = '-f alsa -ac 2 -i hw:0,0'
#audio_options = '-f alsa -ac 2 -i hw:0'

# for the c-media (SIIG) usb device
# XXX may have to change the number after reboots. ug.
audio_options = '-f alsa -ac 1 -i hw:1,0'

#cmd = 'ffmpeg ' + audio_options + ' -f x11grab -s ' + size + ' -r 30 -i :0.0  -target ntsc-dvd ' + outfile
cmd = 'avconv ' + audio_options + ' -f x11grab -s ' + size + ' -r 30 -i :0.0  -target ntsc-dvd ' + outfile

print cmd
os.system(cmd)
