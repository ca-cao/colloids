import ovito
from ovito.io import *
from ovito.vis import *
import sys
from ovito.modifiers import *

fname = str(sys.argv[1])[0:-4]
node = import_file(str(sys.argv[1]), columns = 
          [None ,"Position.X", "Position.Y", "Position.Z"])
node.add_to_scene()
node.source.cell.display.render_cell = False
node.source.particle_properties.position.display.radius = .5
# esto agrega el color
#modifier = ColorCodingModifier(
#            particle_property = "Potential Energy",
#                gradient = ColorCodingModifier.Magma()
#                )
#node.modifiers.append(modifier)

# Let OVITO render an image of the active viewport.
vp = ovito.dataset.viewports.active_vp
vp.type = Viewport.Type.PERSPECTIVE
#vp.camera_pos = (100, 50, 50)
#vp.camera_dir = (-100, -50, -50)
rs = RenderSettings(
    size = (1024,768),
    background_color = (0,0,0)
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
    #background_color = (0, 0, 0, 0)
)
vp.zoom_all()
image = vp.render(rs)



# Save image to disk.
image.save(fname+"top.png")


exit()


# Paint something on top of the rendered image.
#painter = QPainter(image)
#painter.drawText(10, 20, "Hello world!")
#del painter
