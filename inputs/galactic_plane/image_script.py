from aplpy import FITSFigure
from astropy.io import fits
filename = raw_input("File to plot: ")
hdu = fits.open(filename)[1]
fig = FITSFigure(hdu)
fig.show_colorscale(stretch='log', interpolation='none')
#fig.add_colorbar()
#fig.colorbar.set_ticks([0.1e-11, 0.4e-11, 1.2e-11])
#fig.colorbar.set_width(0.1)
#fig.colorbar.set_font('serif')
#fig.colorbar.set_label_properties(size='small')
#fig.hide_ytick_labels()
#fig.hide_yaxis_label()
#fig.hide_xaxis_label()
fig.tick_labels.set_yformat('ddd')
fig.tick_labels.set_xformat('ddd')
fig.set_tick_yspacing(5)
fig.axis_labels.set_font(size='small', family='roman')
fig.set_tick_labels_font('serif')
fig.set_tick_labels_size('small')
fig.save('image.pdf')

