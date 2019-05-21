import ovito
from ovito.io import *
from ovito.vis import *
import sys
from ovito.modifiers import *
import matplotlib
import matplotlib.pyplot as plt
from PyQt5.QtGui import *

# funcion para poner una barra
def render(painter,modifier ,**args):

	# Find the existing HistogramModifier in the pipeline 
	# and get its histogram data.
	for mod in ovito.dataset.selected_node.modifiers:
		if isinstance(mod, HistogramModifier):
			x = mod.histogram[:,0]
			y = mod.histogram[:,1]
			break
	if not 'x' in locals():
		raise RuntimeError('Histogram modifier not found.')
	
	# Get size of rendered viewport image in pixels.
	viewport_width = painter.window().width()
	viewport_height = painter.window().height()
	
	#  Compute plot size in inches (DPI determines label size)
	dpi = 80
	plot_width = 0.5 * viewport_width / dpi
	plot_height = 0.5 * viewport_height / dpi
	
	# Create figure
	fig, ax = plt.subplots(figsize=(plot_width,plot_height), dpi=dpi)
	fig.patch.set_alpha(0.5)
	plt.title('Hexatic order parameter')
	
	# Plot histogram data
	ax.bar(x, y)
	plt.tight_layout()
	
	# Render figure to an in-memory buffer.
	buf = fig.canvas.print_to_buffer()
	
	# Create a QImage from the memory buffer
	res_x, res_y = buf[1]
	img = PyQt5.QtGui.QImage(buf[0], res_x, res_y, PyQt5.QtGui.QImage.Format_RGBA8888)
	
	# Paint QImage onto rendered viewport 
	painter.drawImage(0,0,img)



fname = str(sys.argv[1])[0:-4]
node = import_file(str(sys.argv[1]), columns = 
          ["Position.X", "Position.Y", "Position.Z","Potential Energy"])
node.add_to_scene()
node.source.cell.display.render_cell = False
node.source.particle_properties.position.display.radius = .56
# esto agrega el color
modifier = ColorCodingModifier(
            particle_property = "Potential Energy",
                gradient = ColorCodingModifier.Magma()
                )
node.modifiers.append(modifier)

# Let OVITO render an image of the active viewport.
vp = ovito.dataset.viewports.active_vp
vp.type = Viewport.Type.PERSPECTIVE
#vp.camera_pos = (100, 50, 50)
#vp.camera_dir = (-100, -50, -50)
rs = RenderSettings(
    size = (1024,768),
    #background_color = (0,0,0)
)
vp.zoom_all()
image = vp.render(rs)


# Save image to disk.
image.save(fname+"prsp.png")


vp.type = Viewport.Type.TOP
#vp.camera_pos = (100, 50, 50)
#vp.camera_dir = (-100, -50, -50)
rs = RenderSettings(
    size = (1024,768),
    background_color = (0,0,0)
)
vp.zoom_all()
image = vp.render(rs)


# barra de color
painter = QPainter(image)
render(painter,modifier)

# Save image to disk.
image.save(fname+"top.png")


exit()


# Paint something on top of the rendered image.
#painter = QPainter(image)
#painter.drawText(10, 20, "Hello world!")
#del painter
