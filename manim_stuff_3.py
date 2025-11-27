from manim import *
from manim.opengl import *
import numpy as np
import sympy as sp

class test_polygon(MovingCameraScene):
    def construct(self):
        # Load 2D coordinates (x, y)
        tooth_star_2d = np.loadtxt("28102025_085750.txt")
        peak_star_2d = np.loadtxt("28102025_110553.txt")
        corner_square_2d = np.loadtxt("28102025_111212.txt")

        # Convert to 3D by adding a z=0 coordinate
        tooth_star_3d = np.column_stack([tooth_star_2d, np.zeros(len(tooth_star_2d))])
        peak_star_3d = np.column_stack([peak_star_2d, np.zeros(len(peak_star_2d))])
        corner_square_3d = np.column_stack([corner_square_2d, np.zeros(len(corner_square_2d))])

        # Create the polygon
        tooth_star = Polygon(*tooth_star_3d, fill_opacity=0, color=WHITE)
        peak_star = Polygon(*peak_star_3d, fill_opacity=0, color=WHITE)
        corner_square_tr = Polygon(*corner_square_3d, fill_opacity=0, color=WHITE)
        corner_square_tl = corner_square_tr.copy().rotate_about_origin(angle=90*DEGREES)

        self.camera.frame_height = 40
        self.camera.frame_width = 40
        
        # Show it on screen
        self.play(
            Create(tooth_star),
            Create(peak_star),
            run_time = 2
            )
        

        self.play(
            tooth_star.animate.scale(2).rotate_about_origin(angle=-45*DEGREES),
            peak_star.animate.scale(0.4).rotate_about_origin(angle=90*DEGREES),
            self.camera.frame.animate.set(width=17.5),
            run_time=2
            )
        
        self.play(
            self.camera.frame.animate.set(width=20),
            tooth_star.animate.scale(0.25),
            peak_star.animate.scale(1/0.4).rotate_about_origin(angle=-180*DEGREES),
        )

        combined_stars = VGroup(tooth_star,peak_star)

        self.play(
            Rotate(combined_stars,angle=180*DEGREES),
            self.camera.frame.animate.set(width=10),
            run_time = 2
        )



        self.wait()

class ScatterWithSlice(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75*DEGREES, theta=-45*DEGREES)
        x_data = np.loadtxt("det_x_export.txt")
        y_data = np.loadtxt("det_y_export.txt")
        TOF_data = np.loadtxt("det_TOF_export.txt")

        combined = [x_data[::3],y_data[::3],TOF_data[::3]]
        points3d = np.transpose(combined)
        
        print("Files loaded")


        # --- Create left-side 3D scatter ---
        axes3d = ThreeDAxes(
            x_range=[-0.02, 0.02, 0.005],  # Range from -10 to 10, with steps of 2
            y_range=[-0.02, 0.02, 0.005],   # Range from -5 to 5, with steps of 1
            z_range=[100, 135, 15],    # Range from 0 to 8, with steps of 2
            x_length=10,           # Adjust the visual length of the x-axis
            y_length=10,           # Adjust the visual length of the y-axis
            z_length=4            # Adjust the visual length of the z-axis
        )

        self.add(axes3d)

        for i in np.arange(len(points3d)):

            self.add(Dot3D(point=axes3d.coords_to_point(points3d[i,0],points3d[i,1],1e9*points3d[i,2]),radius=0.01,color=WHITE))



        self.wait(2)

class ScatterWithSlice2(ThreeDScene):
    def construct(self):
        #self.set_camera_orientation(phi=75*DEGREES, theta=-45*DEGREES)

        # Load data exactly as you are doing
        x_data = np.loadtxt("det_x_export.txt")
        y_data = np.loadtxt("det_y_export.txt")
        TOF_data = np.loadtxt("det_TOF_export.txt")

        # Downsample exactly your way
        combined = [x_data[::3], y_data[::3], TOF_data[::3]]
        points3d = np.transpose(combined)

        # --- Create left 3D axes the same way ---
        axes3d = ThreeDAxes(
            x_range=[-0.02, 0.02, 0.005],
            y_range=[-0.02, 0.02, 0.005],
            z_range=[100, 135, 5],
            x_length=10,
            y_length=10,
            z_length=4
        ).shift(LEFT*4).rotate(phi=75*DEGREES, theta=-45*DEGREES)

        self.add(axes3d)

        # --- Plot points EXACTLY the way you currently do ---
        for i in range(len(points3d)):
            self.add(Dot3D(
                point=axes3d.c2p(points3d[i,0], points3d[i,1], 1e9*points3d[i,2]),
                radius=0.01,
                color=WHITE
            ))

        # --- Slice plane ---
        plane_x = ValueTracker(0.0)

        plane = Rectangle(width=0.002, height=10, fill_opacity=0.3, fill_color=YELLOW)
        plane.rotate(PI/2, axis=UP)  # vertical plane
        plane.add_updater(
            lambda m: m.move_to(axes3d.c2p(plane_x.get_value(), 0, 115))
        )
        self.add(plane)

        # --- Right side 2D axes ---
        axes2d = Axes(
            x_range=[-0.02, 0.02, 0.005],
            y_range=[100, 135, 5],
            x_length=6,
            y_length=6
        ).shift(RIGHT*4)

        self.add(axes2d)

        # --- Group where the 2D slice points will go ---
        slice_2d_group = VGroup().shift(RIGHT*4)
        self.add(slice_2d_group)

        SLICE_THICKNESS = 0.002  # thickness in x

        # --- Live updater function ---
        def update_slice(group):
            group.submobjects.clear()
            x_val = plane_x.get_value()

            # Select points near the plane
            mask = np.abs(points3d[:,0] - x_val) < SLICE_THICKNESS
            slice_pts = points3d[mask]

            for px, py, ptof in slice_pts:
                group.add(
                    Dot(
                        axes2d.c2p(py, 1e9*ptof),  # (y,z) plane
                        radius=0.06,
                        color=BLUE
                    )
                )

        slice_2d_group.add_updater(update_slice)

        # --- Animate slicing ---
        self.play(plane_x.animate.set_value(0.02), run_time=3)
        self.play(plane_x.animate.set_value(-0.02), run_time=3)

        self.wait(2)


class TOF_scene(ThreeDScene):
    def construct(self):
        x_data = np.loadtxt("det_x_export.txt")
        y_data = np.loadtxt("det_y_export.txt")
        TOF_data = np.loadtxt("det_TOF_export.txt")

        combined = [x_data[::100],y_data[0::100],TOF_data[0::100]]
        points3d = np.transpose(combined)

        self.set_camera_orientation(phi=75 * DEGREES, theta=-30 * DEGREES)

        TOF_plane = Surface(
            points3d,
            resolution=(24,24),
        )

        self.add(TOF_plane)
        self.wait(2)

class TOF_scene2(ThreeDScene):
    def construct(self):
        x = np.loadtxt("det_x_export.txt")
        y = np.loadtxt("det_y_export.txt")
        z = np.loadtxt("det_TOF_export.txt")

        points3d = np.column_stack([x, y, z])

        self.set_camera_orientation(phi=75*DEGREES, theta=-30*DEGREES)

        axes = ThreeDAxes()
        
        # The magic line: one object, many points
        cloud = OpenGLPointCloud(
            points3d,
            point_radius = 0.03,       # same visual size as Dot3D
            color = WHITE
        )

        self.add(axes, cloud)
        self.wait(2)

class const_line(MovingCameraScene):
    def construct(self):
        self.add_anchor_and_orbit()
        self.add_rotor_and_lines()

    
    # Add the line anchor on the right side
    def add_anchor_and_orbit(self):
        orbit_radius = 5
        self.origin_point = np.array([0,0,0])
        anchor = Dot([orbit_radius,0,0],radius=0.05,color=MAROON)
        self.add(anchor)
        self.anchor = anchor
        #anchor_orbit = Square(side_length=10,color=BLUE)
        anchor_orbit = Polygon(
            [5,5,0],
            [-5,5,0],
            [-5,-5,0],
            [5,-5,0],
            color=BLUE
        )
        self.anchor_orbit = anchor_orbit
        self.add(anchor_orbit)
        anchor.move_to(anchor_orbit.point_from_proportion(0))
        self.t_offset = 0
        rate = -0.3

        def anchor_rotation(mob,dt):
            self.t_offset += (dt * rate)
            mob.move_to(anchor_orbit.point_from_proportion(self.t_offset % 1))

        def square_rotation(mob,dt):
            mob.rotate(dt * -0.5)


        anchor.add_updater(anchor_rotation)
        anchor_orbit.add_updater(square_rotation)

        
        self.collected_points = []
        self.target_polygon = anchor_orbit
        self.ready_to_morph = False   # signal flag

        # Helper to test closeness
        def near(a, b, eps=0.015):
            return np.linalg.norm(a - b) < eps

        # Updater 1: Detect when anchor reaches vertex
        def corner_detector(mob,dt):
            ac = anchor.get_center()
            corners = self.target_polygon.get_vertices()

            for c in corners:
                if near(ac, c):
                    # avoid adding the same vertex many frames in a row
                    if len(self.collected_points) == 0 or not np.allclose(c, self.collected_points[-1]):
                        self.collected_points.append(c)
                        #self.add(Dot(self.collected_points[-1],radius=0.05,color=WHITE))

                    if len(self.collected_points) == 4:
                        self.ready_to_morph = True  # signal morph
                        self.old_polygon_vertices = self.collected_points
                        print("Set self.ready_to_morph = True")
                    return

        anchor.add_updater(corner_detector)

        self.morph_parameter = ValueTracker(0)

        def morph_controller(mob, dt):
            if self.ready_to_morph:
                #self.ready_to_morph = False
                print("Morph controller triggered")
                print("Length of collected points:", len(self.collected_points))
                print("Collected points:", self.collected_points)

                # Convert collected points to NumPy array BEFORE clearing
                new_corners = np.array(self.collected_points)
                old_corners = self.target_polygon.get_vertices()[:4]  # take only 4 main corners

                
                # Interpolation parameter
                t = self.morph_parameter.get_value() + dt / 1.0
                t = min(t, 1)

                # Interpolate positions
                interp_corners = old_corners * (1 - t) + new_corners * t

                # Update polygon
                self.target_polygon.become(Polygon(*interp_corners, color=BLUE))
                self.morph_parameter.set_value(t)

                if t >= 1:
                    self.target_polygon.become(Polygon(*new_corners, color=BLUE))
                    self.morph_parameter.set_value(0)
                    self.ready_to_morph = False
                    self.collected_points = []  # now safe to clear
       
        anchor.add_updater(morph_controller)



    def add_rotor_and_lines(self):

        anchor = self.anchor
        fan_blade_top = Polygon(
            [0,0,0],
            [-0.15,0.2,0],
            [0,2,0],
            [0.15,0.2,0],
            color=WHITE,fill_opacity=0,stroke_width=2
        )
        fan_blade_left = fan_blade_top.copy().rotate(angle=90*DEGREES,about_point=ORIGIN)
        fan_blade_down = fan_blade_top.copy().rotate(angle=180*DEGREES,about_point=ORIGIN)
        fan_blade_right = fan_blade_top.copy().rotate(angle=270*DEGREES,about_point=ORIGIN)



        # Add dots to the peaks of the fan blades to track the motion later
        dot_top = Dot([0,2,0],radius=0.05,color=RED)
        dot_left = Dot([-2,0,0],radius=0.05,color=RED)
        dot_down = Dot([0,-2,0],radius=0.05,color=RED)
        dot_right = Dot([2,0,0],radius=0.05,color=RED)

        # Add the lines
        def top_to_anchor():
            return Line(dot_top.get_center(),anchor.get_center()).set_color(MAROON)
        
        def left_to_anchor():
            return Line(dot_left.get_center(),anchor.get_center()).set_color(MAROON)
        
        def down_to_anchor():
            return Line(dot_down.get_center(),anchor.get_center()).set_color(MAROON)
        
        def right_to_anchor():
            return Line(dot_right.get_center(),anchor.get_center()).set_color(MAROON)


        rotor = VGroup(fan_blade_top,fan_blade_left,fan_blade_down,fan_blade_right,dot_top,dot_left,dot_down,dot_right)


        def rotation_updater_rotor(mob,dt):
            mob.rotate(angle=dt)

        self.t_offset_rotor = 0
        rate_rotor = 0.05

        rotor_orbit = Circle(radius=2)


        def rotor_precession(mob,dt):
            self.t_offset_rotor += (dt * rate_rotor)
            mob.move_to(rotor_orbit.point_from_proportion(self.t_offset_rotor % 1))

        rotor.add_updater(rotation_updater_rotor)
        rotor.add_updater(rotor_precession)
        

        top_to_anchor_line = always_redraw(top_to_anchor)
        left_to_anchor_line = always_redraw(left_to_anchor)
        down_to_anchor_line = always_redraw(down_to_anchor)
        right_to_anchor_line = always_redraw(right_to_anchor)

        self.add(rotor)
        self.add(top_to_anchor_line,left_to_anchor_line,down_to_anchor_line,right_to_anchor_line)
        self.wait(15)

class MoveIndividualVertices(Scene):
    def construct(self):
        # Create a square
        square = Square()
        square_copy = square.copy().rotate(angle=40*DEGREES)
        square_copy_vertices = square_copy.get_vertices()

        square.set_fill(BLUE, opacity=0.5)
        self.play(Create(square))
        self.play(Create(Polygon(*square_copy_vertices)))
        self.wait()


        print(
            square.get_vertices()
        )

        self.wait()

class const_line_backup2(MovingCameraScene):
    def construct(self):
        self.add_anchor_and_orbit()
        self.add_rotor_and_lines()
        

    # Add the line anchor on the right side
    def add_anchor_and_orbit(self):
        orbit_radius = 5
        self.origin_point = np.array([0,0,0])
        anchor = Dot([orbit_radius,0,0],radius=0.05,color=MAROON)
        self.add(anchor)
        self.anchor = anchor
        anchor_orbit = Square(side_length=10,color=BLUE)
        self.anchor_orbit = anchor_orbit
        self.add(anchor_orbit)
        anchor.move_to(anchor_orbit.point_from_proportion(0))
        self.t_offset = 0
        rate = -0.1

        def anchor_rotation(mob,dt):
            self.t_offset += (dt * rate)
            mob.move_to(anchor_orbit.point_from_proportion(self.t_offset % 1))

        def square_rotation(mob,dt):
            mob.rotate(dt * -0.5)


        anchor.add_updater(anchor_rotation)
        anchor_orbit.add_updater(square_rotation)

    def add_rotor_and_lines(self):

        anchor = self.anchor
        fan_blade_top = Polygon(
            [0,0,0],
            [-0.15,0.2,0],
            [0,2,0],
            [0.15,0.2,0],
            color=WHITE,fill_opacity=0,stroke_width=2
        )
        fan_blade_left = fan_blade_top.copy().rotate(angle=90*DEGREES,about_point=ORIGIN)
        fan_blade_down = fan_blade_top.copy().rotate(angle=180*DEGREES,about_point=ORIGIN)
        fan_blade_right = fan_blade_top.copy().rotate(angle=270*DEGREES,about_point=ORIGIN)


        # Add dots to the peaks of the fan blades to track the motion later
        dot_top = Dot([0,2,0],radius=0.05,color=RED)
        dot_left = Dot([-2,0,0],radius=0.05,color=RED)
        dot_down = Dot([0,-2,0],radius=0.05,color=RED)
        dot_right = Dot([2,0,0],radius=0.05,color=RED)

        # Add the lines
        def top_to_anchor():
            return Line(dot_top.get_center(),anchor.get_center()).set_color(MAROON)
        
        def left_to_anchor():
            return Line(dot_left.get_center(),anchor.get_center()).set_color(MAROON)
        
        def down_to_anchor():
            return Line(dot_down.get_center(),anchor.get_center()).set_color(MAROON)
        
        def right_to_anchor():
            return Line(dot_right.get_center(),anchor.get_center()).set_color(MAROON)


        rotor = VGroup(fan_blade_top,fan_blade_left,fan_blade_down,fan_blade_right,dot_top,dot_left,dot_down,dot_right)


        def rotation_updater_rotor(mob,dt):
            mob.rotate(angle=dt)

        self.t_offset_rotor = 0
        rate_rotor = 0.05

        rotor_orbit = Circle(radius=2)


        def rotor_precession(mob,dt):
            self.t_offset_rotor += (dt * rate_rotor)
            mob.move_to(rotor_orbit.point_from_proportion(self.t_offset_rotor % 1))

        rotor.add_updater(rotation_updater_rotor)
        rotor.add_updater(rotor_precession)
        

        top_to_anchor_line = always_redraw(top_to_anchor)
        left_to_anchor_line = always_redraw(left_to_anchor)
        down_to_anchor_line = always_redraw(down_to_anchor)
        right_to_anchor_line = always_redraw(right_to_anchor)

        self.add(rotor)
        self.add(top_to_anchor_line,left_to_anchor_line,down_to_anchor_line,right_to_anchor_line)
        self.wait(60)



class const_line_backup(MovingCameraScene):
    def construct(self):
        self.add_anchor_and_orbit()
        self.add_rotor_and_lines()
        

    # Add the line anchor on the right side
    def add_anchor_and_orbit(self):
        orbit_radius = 5
        self.origin_point = np.array([0,0,0])
        anchor = Dot([orbit_radius,0,0],radius=0.05,color=MAROON)
        self.add(anchor)
        self.anchor = anchor
        anchor_orbit = Circle(radius=orbit_radius)
        self.anchor_orbit = anchor_orbit
        self.add(anchor_orbit)
        anchor.move_to(anchor_orbit.point_from_proportion(0))
        self.t_offset = 0
        rate = -0.1

        def anchor_rotation(mob,dt):
            self.t_offset += (dt * rate)
            mob.move_to(anchor_orbit.point_from_proportion(self.t_offset % 1))

        anchor.add_updater(anchor_rotation)

    def add_rotor_and_lines(self):

        anchor = self.anchor
        fan_blade_top = Polygon(
            [0,0,0],
            [-0.15,0.2,0],
            [0,2,0],
            [0.15,0.2,0],
            color=WHITE,fill_opacity=0,stroke_width=2
        )
        fan_blade_left = fan_blade_top.copy().rotate(angle=90*DEGREES,about_point=ORIGIN)
        fan_blade_down = fan_blade_top.copy().rotate(angle=180*DEGREES,about_point=ORIGIN)
        fan_blade_right = fan_blade_top.copy().rotate(angle=270*DEGREES,about_point=ORIGIN)


        # Add dots to the peaks of the fan blades to track the motion later
        dot_top = Dot([0,2,0],radius=0.05,color=RED)
        dot_left = Dot([-2,0,0],radius=0.05,color=RED)
        dot_down = Dot([0,-2,0],radius=0.05,color=RED)
        dot_right = Dot([2,0,0],radius=0.05,color=RED)

        # Add the lines
        def top_to_anchor():
            return Line(dot_top.get_center(),anchor.get_center()).set_color(MAROON)
        
        def left_to_anchor():
            return Line(dot_left.get_center(),anchor.get_center()).set_color(MAROON)
        
        def down_to_anchor():
            return Line(dot_down.get_center(),anchor.get_center()).set_color(MAROON)
        
        def right_to_anchor():
            return Line(dot_right.get_center(),anchor.get_center()).set_color(MAROON)


        rotor = VGroup(fan_blade_top,fan_blade_left,fan_blade_down,fan_blade_right,dot_top,dot_left,dot_down,dot_right)


        def rotation_updater_rotor(mob,dt):
            mob.rotate(angle=dt)

        self.t_offset_rotor = 0
        rate_rotor = 0.05

        rotor_orbit = Circle(radius=2)


        def rotor_precession(mob,dt):
            self.t_offset_rotor += (dt * rate_rotor)
            mob.move_to(rotor_orbit.point_from_proportion(self.t_offset_rotor % 1))

        rotor.add_updater(rotation_updater_rotor)
        rotor.add_updater(rotor_precession)
        

        top_to_anchor_line = always_redraw(top_to_anchor)
        left_to_anchor_line = always_redraw(left_to_anchor)
        down_to_anchor_line = always_redraw(down_to_anchor)
        right_to_anchor_line = always_redraw(right_to_anchor)

        self.add(rotor)
        self.add(top_to_anchor_line,left_to_anchor_line,down_to_anchor_line,right_to_anchor_line)
        self.wait(60)


class SineCurveUnitCircle(Scene):
    # contributed by heejin_park, https://infograph.tistory.com/230
    def construct(self):
        self.show_axis()
        self.show_circle()
        self.move_dot_and_draw_curve()
        self.wait()

    def show_axis(self):
        x_start = np.array([-6,0,0])
        x_end = np.array([6,0,0])

        y_start = np.array([-4,-2,0])
        y_end = np.array([-4,2,0])

        x_axis = Line(x_start, x_end)
        y_axis = Line(y_start, y_end)

        self.add(x_axis, y_axis)
        self.add_x_labels()

        self.origin_point = np.array([-4,0,0])
        self.curve_start = np.array([-3,0,0])

    def add_x_labels(self):
        x_labels = [
            MathTex(r"\pi"), MathTex(r"2 \pi"),
            MathTex(r"3 \pi"), MathTex(r"4 \pi"),
        ]

        for i in range(len(x_labels)):
            x_labels[i].next_to(np.array([-1 + 2*i, 0, 0]), DOWN)
            self.add(x_labels[i])

    def show_circle(self):
        circle = Circle(radius=1)
        circle.move_to(self.origin_point)
        self.add(circle)
        self.circle = circle

    def move_dot_and_draw_curve(self):
        orbit = self.circle
        origin_point = self.origin_point

        dot = Dot(radius=0.08, color=YELLOW)
        dot.move_to(orbit.point_from_proportion(0))
        self.t_offset = 0
        rate = 0.25

        def go_around_circle(mob, dt):
            self.t_offset += (dt * rate)
            # print(self.t_offset)
            mob.move_to(orbit.point_from_proportion(self.t_offset % 1))

        def get_line_to_circle():
            return Line(origin_point, dot.get_center(), color=BLUE)

        def get_line_to_curve():
            x = self.curve_start[0] + self.t_offset * 4
            y = dot.get_center()[1]
            return Line(dot.get_center(), np.array([x,y,0]), color=YELLOW_A, stroke_width=2 )


        self.curve = VGroup()
        self.curve.add(Line(self.curve_start,self.curve_start))
        def get_curve():
            last_line = self.curve[-1]
            x = self.curve_start[0] + self.t_offset * 4
            y = dot.get_center()[1]
            new_line = Line(last_line.get_end(),np.array([x,y,0]), color=YELLOW_D)
            self.curve.add(new_line)

            return self.curve

        dot.add_updater(go_around_circle)

        origin_to_circle_line = always_redraw(get_line_to_circle)
        dot_to_curve_line = always_redraw(get_line_to_curve)
        sine_curve_line = always_redraw(get_curve)

        self.add(dot)
        self.add(orbit, origin_to_circle_line, dot_to_curve_line, sine_curve_line)
        self.wait(8.5)

        dot.remove_updater(go_around_circle)


class TexTest(Scene):
    def construct(self):
        text = Tex(R"$\frac{1}{2}=\frac{2}{4}$")

        self.add(text)
        self.wait(3)

class sympy_trials(MovingCameraScene):
    def construct(self):
        self.camera.frame_width = 10
        self.camera.frame_height = 10

        # Time tracker
        self.t_offset = 0

        # Event trackers
        self.event_triggered = False

        # Parameter definitions
        square_width_from_center = 3
        square_rotation_rate = 2
        square_growth_rate = 0.1

        anchor_radius = 0.1
        anchor_movement_rate = 0.25


        # Shape definitions
        
        square = Square(side_length=2*square_width_from_center,fill_opacity=0,color=WHITE)
        anchor = Dot(radius=anchor_radius,color=WHITE).move_to(square.point_from_proportion(0))
        anchor_trace = TracedPath(lambda: anchor.get_center(),stroke_width=3,stroke_color='#ff77ff')


        # Updater definitions
        def square_time_tracker(mob,dt):
            self.t_offset += dt

        def square_size_updater(mob,dt):
            mob.scale(1+dt*square_growth_rate)

        def square_rotation_updater(mob, dt):
            mob.rotate(angle=dt*square_rotation_rate)
            #if self.t_offset >= 1 and self.t_offset < 1 + dt:
            if self.t_offset == 1:
                self.add(Square(side_length=1,color=RED).rotate(angle=self.t_offset*5))

        def pause_event_updater(dt):
            if not self.event_triggered and self.t_offset >= 1:
                self.event_triggered = True

        self.add_updater(pause_event_updater)

        square.add_updater(square_time_tracker)
        #square.add_updater(square_size_updater)
        square.add_updater(square_rotation_updater)
        
        def anchor_move_from_proportion_updater(mob,dt):
            mob.move_to(square.point_from_proportion(anchor_movement_rate*self.t_offset % 1))
            
        anchor.add_updater(anchor_move_from_proportion_updater)

        rotary_group_1 = VGroup(square,anchor)


        # Add mobjects to scene
        self.add(rotary_group_1)
        self.add(anchor_trace)


        self.wait(10)




class sympy_trials_chatGPT(MovingCameraScene):
    def construct(self):
        self.camera.frame_width = 10
        self.camera.frame_height = 10

        # Time tracker
        self.t_offset = 0
        self.event_triggered = False  # flag so we only trigger once

        # Parameters
        square_width_from_center = 3
        square_rotation_rate = 2
        anchor_radius = 0.1
        anchor_movement_rate = 0.25

        # Shapes
        square = Square(side_length=2*square_width_from_center, fill_opacity=0, color=WHITE)
        anchor = Dot(radius=anchor_radius, color=WHITE).move_to(square.point_from_proportion(0))
        anchor_trace = TracedPath(lambda: anchor.get_center(), stroke_width=3, stroke_color='#ff77ff')

        # === Updaters ===

        # Time tracker
        def square_time_tracker(mob, dt):
            self.t_offset += dt

        # Rotation
        def square_rotation_updater(mob, dt):
            mob.rotate(angle=dt * square_rotation_rate)

        # Anchor motion
        def anchor_move_from_proportion_updater(mob, dt):
            mob.move_to(square.point_from_proportion(anchor_movement_rate * self.t_offset % 1))

        # Controller for timed event
        def pause_event_updater(dt):
            if not self.event_triggered and self.t_offset >= 1.0:
                self.event_triggered = True  # ensure runs only once

                # Pause movement
                square.remove_updater(square_rotation_updater)
                anchor.remove_updater(anchor_move_from_proportion_updater)

                # Run your event animation
                self.play(Create(Square(side_length=1, color=RED, fill_opacity=0)))

                # Resume movement
                square.add_updater(square_rotation_updater)
                anchor.add_updater(anchor_move_from_proportion_updater)

        # Attach updaters
        square.add_updater(square_time_tracker)
        square.add_updater(square_rotation_updater)
        anchor.add_updater(anchor_move_from_proportion_updater)
        self.add_updater(pause_event_updater)

        # Add objects
        self.add(square, anchor, anchor_trace)

        self.wait(10)

class DFT_rotor(MovingCameraScene):
    def construct(self):
        self.camera.frame_width = 10
        self.camera.frame_height = 10

        # Initialise time tracker variable
        self.t_offset = 0

        # Parameters
        self.blade_height = 2
        self.shrink_factor = 0.4
        self.rotation_rate_1 = 0.5
        self.rotation_rate_2 = 0.7

        # Scene function calls
        self.time_tracker()
        self.add_fan_groups()

        

    # Time tracker
    def time_tracker(self):
        def scene_time_tracker_updater(dt):
            self.t_offset += dt
        self.add_updater(scene_time_tracker_updater)


    def add_fan_groups(self):   
        tracker_dot = Dot([0,0,0],radius=0.05,color=MAROON) 
        fan_blade_1 = Polygon(
            [0,0,0],
            [-0.15,0.2,0],
            [0,self.blade_height,0],
            [0.15,0.2,0],
            color=WHITE,fill_opacity=0,stroke_width=2
        )

        # Add dots to the peaks of the fan blades to track the motion later
        dot_top_1 = Dot([0,self.blade_height,0],radius=0,color=RED)
        self.dot_top_1 = dot_top_1
        dot_bottom_1 = Dot([0,0,0],radius=0,color=RED)
        self.dot_bottom_1 = dot_bottom_1

        fan_group_1 = VGroup(fan_blade_1,dot_top_1,dot_bottom_1)
        fan_group_2 = fan_group_1.copy().scale(self.shrink_factor)
        

        self.fan_group_1 = fan_group_1
        self.fan_group_2 = fan_group_2

        # Align base of small fan with tip of large fan
        base_shift = dot_top_1.get_center() - fan_group_2[2].get_center()
        fan_group_2.shift(base_shift)

        def universal_updater(mob,dt):
            self.fan_group_1.rotate(angle=dt*self.rotation_rate_1,about_point=ORIGIN)
            #self.fan_group_2.move_to(self.dot_top_1.get_center()- self.fan_group_2[2].get_center())
            #self.fan_group_2.rotate(angle=dt*self.rotation_rate_2,about_point=self.dot_top_1.get_center())
            

        tracker_dot.add_updater(universal_updater)

        # Add mobjects
        self.add(fan_group_1,fan_group_2)

        # Export timer
        self.wait(5)

class DFT_rotor_2(MovingCameraScene):
    def construct(self):
        self.camera.frame_width = 8
        self.camera.frame_height = 8

        # Initialise time tracker variable
        self.t_offset = 0

        # Parameters
        self.blade_height = 2
        self.shrink_factor = 0.6
        self.rotation_rate_1 = 1.5
        self.rotation_rate_2 = 4
        self.fan_group_2_shrink_rate = 0.02
        self.group_translation_rate = 5

        # Scene function calls
        self.time_tracker()

        self.fan_blade_1 = Polygon(
            [0,0,0],
            [-0.15,0.2,0],
            [0,self.blade_height,0],
            [0.15,0.2,0],
            color=WHITE,fill_opacity=0.2,stroke_width=2
        )

        # Add dots to the peaks of the fan blades to track the motion later
        self.dot_top_1 = Dot([0,self.blade_height,0],radius=0.05,color=RED)
        self.dot_top_1_coordinate_store = self.dot_top_1.get_center()
        self.dot_bottom_1 = Dot([0,0,0],radius=0.05,color=RED)

        self.fan_group_1 = VGroup(self.fan_blade_1,self.dot_top_1,self.dot_bottom_1)
        self.fan_group_2 = self.fan_group_1.copy().scale(self.shrink_factor)
        
        
        # Align base of small fan with tip of large fan
        base_shift = self.dot_top_1.get_center() - self.fan_group_2[2].get_center()
        self.fan_group_2.shift(base_shift)

        # Updaters
        def fan_group_1_rotation_updater(mob,dt):
            #mob.shift(LEFT*dt*self.group_translation_rate)
            #self.camera.frame.move_to(mob[2].get_center())
            #mob.rotate(angle=dt*self.rotation_rate_1,about_point=ORIGIN)
            mob.rotate(angle=dt*self.rotation_rate_1,about_point=mob[2].get_center())
            self.dot_top_1_coordinate_store = mob[1].get_center()
            

        self.fan_group_1.add_updater(fan_group_1_rotation_updater)

        def fan_group_2_rotation_updater(mob,dt):
            mob.shift(self.dot_top_1_coordinate_store-mob[2].get_center())
            mob.rotate(angle=dt*self.rotation_rate_2,about_point=self.dot_top_1.get_center())

        self.fan_group_2.add_updater(fan_group_2_rotation_updater)

        def fan_group_2_size_updater(mob,dt):
            mob.scale(1+self.fan_group_2_shrink_rate*np.sin(self.t_offset*2*np.pi))

        #self.fan_group_2.add_updater(fan_group_2_size_updater)
        self.pointer_trace = TracedPath(lambda: self.fan_group_2[1].get_center(),stroke_color=[WHITE, RED],stroke_width=2)


        self.add(self.fan_group_1,self.fan_group_2)
        self.add(self.pointer_trace)
        self.wait(60)
        
class DFT_rotor_3(MovingCameraScene):
    def construct(self):
        self.camera.frame_width = 8
        self.camera.frame_height = 8

        # Initialise time tracker variable
        self.t_offset = 0

        # Parameters
        self.blade_height = 2
        self.shrink_factor = 0.6
        self.rotation_rate_1 = 1.5
        self.rotation_rate_2 = 4
        self.fan_group_2_shrink_rate = 0.02
        self.group_translation_rate = 5

        # Scene function calls
        self.time_tracker()


        # Shape definitions
        self.fan_blade_1 = Polygon(
            [0,0,0],
            [-0.15,0.2,0],
            [0,self.blade_height,0],
            [0.15,0.2,0],
            color=WHITE,fill_opacity=0.2,stroke_width=2
        )

        self.orbit_route = Square(side_length=2,color=BLUE,fill_opacity=0)
        self.add(self.orbit_route)

        # Add dots to the peaks of the fan blades to track the motion later
        self.dot_top_1 = Dot([0,self.blade_height,0],radius=0.05,color=RED)
        self.dot_top_1_coordinate_store = self.dot_top_1.get_center()
        self.dot_bottom_1 = Dot([0,0,0],radius=0.05,color=RED)

        self.fan_group_1 = VGroup(self.fan_blade_1,self.dot_top_1,self.dot_bottom_1)
        self.fan_group_2 = self.fan_group_1.copy().scale(self.shrink_factor)
        
        
        # Align base of small fan with tip of large fan
        base_shift = self.dot_top_1.get_center() - self.fan_group_2[2].get_center()
        self.fan_group_2.shift(base_shift)

        # Updaters
        def fan_group_1_rotation_updater(mob,dt):
            mob.shift(LEFT*dt*self.group_translation_rate)
            #self.camera.frame.move_to(mob[2].get_center())
            #mob.rotate(angle=dt*self.rotation_rate_1,about_point=ORIGIN)
            mob.rotate(angle=dt*self.rotation_rate_1,about_point=mob[2].get_center())
            self.dot_top_1_coordinate_store = mob[1].get_center()
            

        self.fan_group_1.add_updater(fan_group_1_rotation_updater)

        def fan_group_2_rotation_updater(mob,dt):
            mob.shift(self.dot_top_1_coordinate_store-mob[2].get_center())
            mob.rotate(angle=dt*self.rotation_rate_2,about_point=self.dot_top_1.get_center())

        self.fan_group_2.add_updater(fan_group_2_rotation_updater)

        def fan_group_2_size_updater(mob,dt):
            mob.scale(1+self.fan_group_2_shrink_rate*np.sin(self.t_offset*2*np.pi))

        #self.fan_group_2.add_updater(fan_group_2_size_updater)
        self.pointer_trace = TracedPath(lambda: self.fan_group_2[1].get_center(),stroke_color=[WHITE, RED],stroke_width=2)


        self.add(self.fan_group_1,self.fan_group_2)
        self.add(self.pointer_trace)
        self.wait(60)

    # Time tracker
    def time_tracker(self):
        def scene_time_tracker_updater(dt):
            self.t_offset += dt
        self.add_updater(scene_time_tracker_updater)

class checkerboard(MovingCameraScene):
    def construct(self):
        self.frame_rate = 165
        self.camera.frame_width = 15
        self.camera.frame_height = 15

        # Parameters
        track_point_radius = 0.1
        self.group_rotations_per_second = 0.25
        self.g1_speed = 0.25
        self.g2_speed = 0.25
        self.g3_speed = 0.25
        self.g4_speed = 0.25
        

        # Initialise time tracker variable
        self.t_offset = 0
        self.time_tracker()

        self.square = Square(side_length=2,color=WHITE,fill_opacity=0)
        self.dot = Dot(color=MAROON,radius=track_point_radius).move_to(self.square.point_from_proportion(0))
        self.group = VGroup(self.square,self.dot)

        self.group1 = self.group.copy().shift(UP*2+LEFT*2)
        self.group2 = self.group.copy().shift(UP*2+RIGHT*2)
        self.group3 = self.group.copy().shift(DOWN*2+RIGHT*2)
        self.group4 = self.group.copy().shift(DOWN*2+LEFT*2)

        #globals().update(locals())
        def group_updater_1(mob,dt):
            mob[0].rotate(angle = self.group_rotations_per_second*dt * 2 * np.pi)
            #mob[1].move_to(mob[0].point_from_proportion(self.g1_speed*self.t_offset % 1))
            mob[1].move_to(mob[0].point_from_proportion(0))

        def group_updater_2(mob,dt):
            mob[0].rotate(angle = self.group_rotations_per_second*dt * 2 * np.pi)
            #mob[1].move_to(mob[0].point_from_proportion(self.g2_speed*self.t_offset % 1))
            mob[1].move_to(mob[0].point_from_proportion(0.75))

        def group_updater_3(mob,dt):
            mob[0].rotate(angle = self.group_rotations_per_second*dt * 2 * np.pi)
            #mob[1].move_to(mob[0].point_from_proportion(self.g3_speed*self.t_offset % 1))
            mob[1].move_to(mob[0].point_from_proportion(0.5))

        def group_updater_4(mob,dt):
            mob[0].rotate(angle = self.group_rotations_per_second*dt * 2 * np.pi)
            #mob[1].move_to(mob[0].point_from_proportion(self.g4_speed*self.t_offset % 1))
            mob[1].move_to(mob[0].point_from_proportion(0.25))

            

        self.group1.add_updater(group_updater_1)
        self.group2.add_updater(group_updater_2)
        self.group3.add_updater(group_updater_3)
        self.group4.add_updater(group_updater_4)
        self.add(self.group1,self.group2,self.group3,self.group4)

        def line_group1_group2():
            return Line(self.group1[1].get_center(), self.group2[1].get_center(), color=MAROON)
        
        def line_group2_group3():
            return Line(self.group2[1].get_center(), self.group3[1].get_center(), color=MAROON)
        
        def line_group3_group4():
            return Line(self.group3[1].get_center(), self.group4[1].get_center(), color=MAROON)
        
        def line_group4_group1():
            return Line(self.group4[1].get_center(), self.group1[1].get_center(), color=MAROON)
        
        line1_2 = always_redraw(line_group1_group2)
        line2_3 = always_redraw(line_group2_group3)
        line3_4 = always_redraw(line_group3_group4)
        line4_1 = always_redraw(line_group4_group1)

        self.add(line1_2,line2_3,line3_4,line4_1)


        self.wait(30)



    # Time tracker
    def time_tracker(self):
        def scene_time_tracker_updater(dt):
            self.t_offset += dt
        self.add_updater(scene_time_tracker_updater)


class spectral_phase(MovingCameraScene):
    def construct(self):
        aspect_ratio_scale = 2
        self.camera.frame_width = 19.20*aspect_ratio_scale
        self.camera.frame_height = 10.80*aspect_ratio_scale

        # Add time tracker to scene
        self.t_offset = 0
        #self.time_tracker()

        # Parameters
        x_start = -25
        x_end = 25
        frequency_global_scale = 10
        x01_start_value = -14
        x02_start_value = -6
        x03_start_value = -2 
        CEP_1_start_value = 0
        CEP_2_start_value = np.pi/4
        CEP_3_start_value = np.pi/6


        frequency_1 = 0.8 # oscillations per unit distance in x
        amplitude_1 = 1
        CEP_1 = ValueTracker(CEP_1_start_value)
        x0_1 = ValueTracker(x01_start_value)
        sigma_1 = 0.8

        frequency_2 = 1.3 # oscillations per unit distance in x
        amplitude_2 = 1.5
        CEP_2 = ValueTracker(CEP_2_start_value)
        x0_2 = ValueTracker(x02_start_value)
        sigma_2 = 1.2

        frequency_3 = 2 # oscillations per unit distance in x
        amplitude_3 = 1.2
        CEP_3 = ValueTracker(CEP_3_start_value)
        x0_3 = ValueTracker(x03_start_value)
        sigma_3 = 1.35

        # Add background number plane
        numberplane = NumberPlane(
            x_range = [-50,50],
            y_range = [-10,10],
            x_length = 100,
            y_length= 20
        )

        self.add(numberplane)

        # Function that will give envelopes
        def gaussian(x, A=1, x0=0, sigma=1):
            return A * np.exp(-(x - x0)**2 / (2 * sigma**2))
        
        def get_field_graph(x, frequency, A, x0, sigma):
            return np.cos(frequency * x * 2 * np.pi)*gaussian(x, A, x0, sigma)
        
        field_1_graph = FunctionGraph(
            lambda x: get_field_graph((x+CEP_1.get_value()), frequency_1, amplitude_1, x0_1.get_value(), sigma_1),
            x_range=[x_start, x_end],
            color=BLUE
        )

        field_2_graph = FunctionGraph(
            lambda x: get_field_graph((x+CEP_2.get_value()), frequency_2, amplitude_2, x0_2.get_value(), sigma_2),
            x_range=[x_start, x_end],
            color=MAROON
        )

        field_3_graph = FunctionGraph(
            lambda x: get_field_graph((x+CEP_3.get_value()), frequency_3, amplitude_3, x0_3.get_value(), sigma_3),
            x_range=[x_start, x_end],
            color=GREEN
        )

        # Determine velocities of field graphs based on their frequency
        fg1_velocity = frequency_global_scale*1/frequency_1
        fg2_velocity = frequency_global_scale*1/frequency_2
        fg3_velocity = frequency_global_scale*1/frequency_3

        def fg1_updater(mob,dt):
            mob.become(
                FunctionGraph(
                    lambda x: get_field_graph((x+CEP_1.get_value()), frequency_1, amplitude_1, x0_1.get_value(), sigma_1),
                    x_range=[x_start, x_end],
                    color=BLUE
                ).shift(UP*3)
            )

        field_1_graph.add_updater(fg1_updater)

        def fg2_updater(mob,dt):
            mob.become(
                FunctionGraph(
                    lambda x: get_field_graph((x+CEP_2.get_value()), frequency_2, amplitude_2, x0_2.get_value(), sigma_2),
                    x_range=[x_start, x_end],
                    color=MAROON
                )
            )

        field_2_graph.add_updater(fg2_updater)

        def fg3_updater(mob,dt):
            mob.become(
                FunctionGraph(
                    lambda x: get_field_graph((x+CEP_3.get_value()), frequency_3, amplitude_3, x0_3.get_value(), sigma_3),
                    x_range=[x_start, x_end],
                    color=GREEN
                ).shift(DOWN*3)
            )

        field_3_graph.add_updater(fg3_updater)



        # Move field graphs so they dont overlap
        field_1_graph.shift(UP*3)
        field_2_graph.shift(UP*0)
        field_3_graph.shift(DOWN*3)

        #self.add(field_1_graph,field_2_graph,field_3_graph)
        
        self.play(
            Create(field_1_graph),
            Create(field_2_graph),
            Create(field_3_graph),
            run_time = 3
        )

        
        # Move standalone graphs in x0 offset to overlap point and back to preview behavior of combined plot
        self.play(
            x0_1.animate.set_value(15),
            x0_2.animate.set_value(15),
            x0_3.animate.set_value(15),
            run_time = 3,
            rate_func = linear
        )
        self.play(
            x0_1.animate.set_value(x01_start_value),
            x0_2.animate.set_value(x02_start_value),
            x0_3.animate.set_value(x03_start_value),
            run_time = 3,
            rate_func = linear
        )
        

        # Create overlapped graph from all field_x_graph components
        field_combined_graph = FunctionGraph(
            lambda x: get_field_graph((x+CEP_1.get_value()), frequency_1, amplitude_1, x0_1.get_value(), sigma_1) + get_field_graph((x+CEP_2.get_value()), frequency_2, amplitude_2, x0_2.get_value(), sigma_2) + get_field_graph((x+CEP_3.get_value()), frequency_3, amplitude_3, x0_3.get_value(), sigma_3),
            x_range=[x_start, x_end],
            color=WHITE
        ) 

        # Create updater for field_combined_graph
        def field_combined_graph_updater(mob,dt):
            mob.become(
                FunctionGraph(
                    lambda x: get_field_graph((x+CEP_1.get_value()), frequency_1, amplitude_1, x0_1.get_value(), sigma_1) + get_field_graph((x+CEP_2.get_value()), frequency_2, amplitude_2, x0_2.get_value(), sigma_2) + get_field_graph((x+CEP_3.get_value()), frequency_3, amplitude_3, x0_3.get_value(), sigma_3),
                    x_range=[x_start, x_end],
                    color=WHITE
                )
            )

        

        for mob in self.mobjects:
            #mob.remove_updater(mob.get_updaters())
            mob.suspend_updating()

        # Play the transform
        field_graph_group = VGroup(field_2_graph,field_1_graph,field_3_graph)
        self.play(
            Transform(
                field_graph_group,field_combined_graph,
                replace_mobject_with_target_in_scene=True
            ),
            run_time = 2
        )

        # Move camera to updater attach position smoothly
        self.play(
            self.camera.frame.animate.move_to([x0_2.get_value(), 0, 0]),
            run_time = 1.5
        )

        # Move camera to valuetracker value of field_grpah_2
        def camera_updater(mob,dt):
            self.camera.frame.move_to([x0_2.get_value(), 0, 0])

        field_combined_graph.add_updater(camera_updater)

        # Add updater to the field combined graph after the transform has replaced the Group of standalone graphs
        field_combined_graph.add_updater(field_combined_graph_updater)

        # Play ValueTrackers again but this time for the combined graph
        self.play(
            x0_1.animate.set_value(15),
            x0_2.animate.set_value(15),
            x0_3.animate.set_value(15),
            run_time = 4,
            rate_func = smooth
        )

        self.play(
            CEP_1.animate.set_value(2*np.pi),
            run_time = 3
        )


        self.wait(2)
    
    # Time tracker
    def time_tracker(self):
        def scene_time_tracker_updater(dt):
            self.t_offset += dt
        self.add_updater(scene_time_tracker_updater)

#config.renderer = "opengl"
#config["preview"] = True   # keep window open after rendering

class opengl_test(MovingCameraScene):
    def construct(self):
        self.t_offset = 0
        self.time_tracker()
        
        
        # Parameters
        square_rotations_per_second = 1

        square = Square(side_length=1,color=WHITE,fill_opacity=0)

        def square_updater(mob,dt):
            mob.rotate(angle=dt*2*np.pi*square_rotations_per_second)

        square.add_updater(square_updater)

        self.add(square)
        self.wait(15)

        self.interactive_embed()

    # Time tracker
    def time_tracker(self):
        def scene_time_tracker_updater(dt):
            self.t_offset += dt
        self.add_updater(scene_time_tracker_updater)


class SwingRotor(MovingCameraScene):
    def construct(self):
        self.camera.frame_width = 30
        self.camera.frame_height = 30

        self.t_offset = 0

        def time_tracker(dt):
            self.t_offset += dt
        self.add_updater(time_tracker)

        # Parameters
        initial_rotor_offset = LEFT*7
        framerate = 165
        blade_height = 2
        rotor_rotation_velocity_scale_factor = 0.25
        target_dot_offset = RIGHT*0
        fan_whole_size_scale_factor = 0.5

        # Shape definitions
        fan_blade_1 = Polygon(
            [0,0,0],
            [-0.15,0.2,0],
            [0,blade_height,0],
            [0.15,0.2,0],
            color=WHITE,fill_opacity=0.2,stroke_width=2
        )

        fan_blade_2  = fan_blade_1.copy().rotate(angle=90*DEGREES,about_point=ORIGIN)
        fan_blade_3  = fan_blade_1.copy().rotate(angle=180*DEGREES,about_point=ORIGIN)
        fan_blade_4  = fan_blade_1.copy().rotate(angle=270*DEGREES,about_point=ORIGIN)

        dot_top_1 = Dot(fan_blade_1.get_vertices()[2],color=RED,radius=0.05)
        dot_top_2 = Dot(fan_blade_2.get_vertices()[2],color=RED,radius=0.05)
        dot_top_3 = Dot(fan_blade_3.get_vertices()[2],color=RED,radius=0.05)
        dot_top_4 = Dot(fan_blade_4.get_vertices()[2],color=RED,radius=0.05)

        fan_group_1 = VGroup(fan_blade_1,dot_top_1)
        fan_group_2 = VGroup(fan_blade_2,dot_top_2)
        fan_group_3 = VGroup(fan_blade_3,dot_top_3)
        fan_group_4 = VGroup(fan_blade_4,dot_top_4)

        fan_whole = VGroup(fan_group_1,fan_group_2,fan_group_3,fan_group_4).scale(fan_whole_size_scale_factor)

        durations = np.array([
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1
        ])
        t_total = sum(durations)

        pts = np.array([
            [1, 1, 0],
            [1, 0, 0],
            [1, 1, 0],
            [3, 1, 0],
            [3, 3, 0],
            [5, 3, 0],
            [5, -2, 0],
            [7, -2, 0],
            [7, 0, 0],
            [8, 0, 0]
        ])

        # Create handles for the cubic bezier path of the target dot
        h1_list, h2_list = get_smooth_cubic_bezier_handle_points(pts)

        # Build all Bézier segments
        beziers = VGroup()
        for i in range(len(pts) - 1):
            b = CubicBezier(
                pts[i],        # anchor start
                h1_list[i],    # handle 1 for segment i
                h2_list[i],    # handle 2 for segment i
                pts[i+1],
                color=RED       # anchor end
            )
            beziers.add(b)
            #self.add(b)

        bezier_sample_range = np.linspace(0,1,1000)
        bezier_function = bezier(pts)
        bezier_path = VMobject()
        bezier_points = bezier_function(bezier_sample_range)
        bezier_path.set_points_smoothly(bezier_points)

        self.add(bezier_path)

        # Numberplane for background
        numberplane = NumberPlane(
            x_range=[-50,50],
            y_range=[-50,50],
            x_length=100,
            y_length=100,
        )

        # Keep numberplane in this configuration & location even if camera angles are changed later
        #last_cam_center = self.camera.frame.get_center()
        #print(last_cam_center)
        def freeze_plane(mob, dt):
            last_cam_center = [0,0,0]
            mob.move_to(last_cam_center)
            #last_cam_center = current

        numberplane.add_updater(freeze_plane)

        # Target for lines from peaks of fanblades
        target_dot = Dot(fan_whole.get_center()+target_dot_offset,color=GREEN,radius=0.1)

        t_tracker = ValueTracker(0)

        def target_dot_updater(mob, dt):
            #alpha = t_tracker.get_value() / t_total
            alpha = min(self.t_offset / t_total, 1)
            mob.move_to(bezier_path.point_from_proportion(alpha))



        target_dot.add_updater(target_dot_updater)
        

        # Add lines from peaks to target
        def top_1_to_target():
            return Line(fan_group_1[1].get_center(),target_dot.get_center(),color=MAROON)
        
        def top_2_to_target():
            return Line(fan_group_2[1].get_center(),target_dot.get_center(),color=MAROON)
        
        def top_3_to_target():
            return Line(fan_group_3[1].get_center(),target_dot.get_center(),color=MAROON)

        def top_4_to_target():
            return Line(fan_group_4[1].get_center(),target_dot.get_center(),color=MAROON)
        
        def target_dot_redraw():
            return target_dot
        
        dot_update = always_redraw(target_dot_redraw)
        line1 = always_redraw(top_1_to_target)
        line2 = always_redraw(top_2_to_target)
        line3 = always_redraw(top_3_to_target)
        line4 = always_redraw(top_4_to_target)

        # Make a wrapper for the rotor rotation that also allows specifying for how long to apply the updater
        def rotor_clockwise_apply(duration,direction):
            #velocity = direction / duration
            def rotor_clockwise(mob,dt):
                mob.elapsed = getattr(mob, "elapsed", 0) + dt
                
                # expose these to the target-dot updater
                mob.current_duration = duration
                mob.current_direction = direction
                if mob.elapsed <= duration:
                    mob.rotate(angle=-dt*2*np.pi/duration * rotor_rotation_velocity_scale_factor)
                    mob.shift(direction/duration * dt)
                    #mob.move_to
                    self.camera.frame.move_to(mob.get_center()-initial_rotor_offset)
                else:
                    mob.remove_updater(rotor_clockwise)
                    del mob.elapsed
            return rotor_clockwise


        self.add(target_dot,numberplane)
        self.add(line1,line2,line3,line4,dot_update)

        # Add shapes

        # Reference dot
        self.add(Dot([1,1,1],radius=0.2,color=RED))

        self.add(fan_whole)

        

        #segment_start_proportion = np.insert(np.cumsum(segment_lengths[:-1]), 0, 0) / total_length
        #segment_end_proportion = np.cumsum(segment_lengths) / total_length
        

        # Movement vectors = differences between consecutive points
        directions = [pts[i+1] - pts[i] for i in range(len(pts)-1)]

        for i in range(len(durations)):
            current_duration = durations[i]
            fan_whole.add_updater(rotor_clockwise_apply(current_duration,directions[i]))
            self.wait(current_duration)

        #self.interactive_embed()
        self.wait(2)

class bezier_test(MovingCameraScene):
    def construct(self):
        self.camera.frame_height = 20
        self.camera.frame_width = 20

        dot1 = Dot([1,1,0],color=WHITE,radius=0.1)
        dot2 = Dot([2,4,0],color=WHITE,radius=0.1)

        self.add(dot1,dot2)

        bezier = CubicBezier(dot1.get_center(),*get_smooth_cubic_bezier_handle_points(dot1.get_center(),dot2.get_center()),dot2.get_center())
        self.play(
            Create(bezier,run_time=3)
        )
        self.wait(1)

class BezierFivePoints(MovingCameraScene):
    def construct(self):
        self.camera.frame_height = 20
        self.camera.frame_width = 20

        # Define 5 points in 2D (z=0)
        pts = [
            np.array([-4, -2, 0]),
            np.array([-1,  2, 0]),
            np.array([ 1, -1, 0]),
            np.array([ 3,  3, 0]),
            np.array([ 5,  0, 0]),
        ]

        # Create dots for visualization
        dots = VGroup(*[Dot(p, color=YELLOW, radius=0.12) for p in pts])
        self.add(dots)

        

        # Compute all cubic handle points  
        # Returns: (handles1, handles2)
        h1_list, h2_list = get_smooth_cubic_bezier_handle_points(pts)

        # Build all Bézier segments
        beziers = VGroup()
        for i in range(len(pts) - 1):
            b = CubicBezier(
                pts[i],        # anchor start
                h1_list[i],    # handle 1 for segment i
                h2_list[i],    # handle 2 for segment i
                pts[i+1]     # anchor end
            )
            beziers.add(b)

        #all_points = np.vstack([curve.get_points() for curve in beziers])
        #bezier_path = VMobject()
        #bezier_path.set_points_smoothly(all_points)

        bezier_sample_range = np.linspace(0,1,1000)
        bezier_function = bezier(pts)
        bezier_path = VMobject()
        bezier_points = bezier_function(bezier_sample_range)
        bezier_path.set_points_smoothly(bezier_points)

        self.add(bezier_path)


        tracking_dot = Dot(radius=0.1,color=RED)
        alpha = ValueTracker(0)

        def tracking_dot_updater(mob,dt):
            mob.move_to(bezier_path.point_from_proportion(alpha.get_value()))

        tracking_dot.add_updater(tracking_dot_updater)
        
        
        

        # Animate drawing all segments
        for curve in beziers:
            self.play(Create(curve, run_time=1.5))
        self.add(tracking_dot)
        self.play(
            alpha.animate.set_value(1),
            run_time = 3
        )

        self.wait(2)

class sublime_test(Scene):
    def construct(self):
        self.add(Circle(radius=1))
        # Checkpoint
        self.add(Square(side_length=8))
        self.add(Square(side_length=7))
        self.add(Circle(radius=2))

        # Checkpoint 2
        self.add(Square(side_length=5,color=GREEN))

        #self.interactive_embed()

class shapeshift_old(MovingCameraScene):
    def construct(self):

        # Parameters
        orbit_rotations_per_second = 0.75

        # Value Trackers
        orbit_angle = ValueTracker(0)
        target_position = ValueTracker(0)

        # Shape definitions
        orbit = Square(side_length=4,color=BLUE,fill_opacity=0)
        orbit_ref = Square(side_length=4,color=BLUE,fill_opacity=0)
        target_dot = Dot(color=MAROON,radius=0.07).move_to(orbit.point_from_proportion(target_position.get_value()))

        # Updaters to change shape based on valuetracker value
        def orbit_updater_wrapper(object):
            reference = object.copy()
            def orbit_updater(mob,dt):
                mob.become(reference.copy().rotate(angle=orbit_angle.get_value()*2*np.pi*orbit_rotations_per_second))

            return orbit_updater

        orbit.add_updater(orbit_updater_wrapper(orbit))

        def target_updater(mob,dt):
            mob.move_to(orbit.point_from_proportion(target_position.get_value()))

        target_dot.add_updater(target_updater)

        self.add(orbit,target_dot)

        new_coords = []
        new_coord_1 = target_dot.copy().move_to(orbit.point_from_proportion(0.25)).rotate(angle=0.25*2*np.pi*orbit_rotations_per_second,about_point=ORIGIN).get_center()
        new_coords.append(new_coord_1)
        new_coord_2 = target_dot.copy().move_to(orbit.point_from_proportion(0.5)).rotate(angle=0.5*2*np.pi*orbit_rotations_per_second,about_point=ORIGIN).get_center()
        new_coords.append(new_coord_2)
        new_coord_3 = target_dot.copy().move_to(orbit.point_from_proportion(0.75)).rotate(angle=0.75*2*np.pi*orbit_rotations_per_second,about_point=ORIGIN).get_center()
        new_coords.append(new_coord_3)
        new_coord_4 = target_dot.copy().move_to(orbit.point_from_proportion(1)).rotate(angle=1*2*np.pi*orbit_rotations_per_second,about_point=ORIGIN).get_center()
        new_coords.append(new_coord_4)

        total_length = 0
        individual_lengths = []
        
        n = len(new_coords)
        for i in range(n):
            segment_length = np.linalg.norm(new_coords[(i+1)%n] - new_coords[i])
            individual_lengths.append(segment_length)
            total_length += segment_length

        proportions = []
        for i in range(len(new_coords)):
            proportions.append(individual_lengths[i]/total_length)

        cumulative_proportions = np.cumsum(proportions)

        print(cumulative_proportions)

        self.add(*[Dot(radius=0.1,color=RED).move_to(newcoord) for newcoord in new_coords])

        self.play(
            orbit_angle.animate.set_value(1),
            target_position.animate.set_value(1),
            rate_func=linear,
            run_time=8
        )
        new_polygon = Polygon(*[coord for coord in new_coords])

        self.play(
            Transform(orbit,new_polygon)
        )

        orbit.remove_updater(orbit.get_updaters())
        orbit.become(new_polygon)
        orbit.add_updater(orbit_updater_wrapper(orbit))

        orbit_angle.set_value(0)
        target_position.set_value(0)

        self.play(
            orbit_angle.animate.set_value(1),
            target_position.animate.set_value(1),
            rate_func=linear,
            run_time=8
        )

class shapeshift(MovingCameraScene):
    def construct(self):
        self.camera.frame_height = 8
        self.camera.frame_width = 8
        # --- Parameters ---
        orbit_rotations_per_second = 1.4
        generation_runtime = 5
        transform_runtime = 1
        new_dot_indication_time = 0.25
        new_dot_radius = 0.08
        tail_dissipation_time = generation_runtime + transform_runtime + new_dot_indication_time
        polygon_stroke_width = 2
        n_generations = 4  # how many iterative polygons to generate
        n_vertices = 4     # number of new vertices per generation. Only tested for "4" so far, good luck

        #Value Trackers
        orbit_angle = ValueTracker(0)
        target_position = ValueTracker(0)

        #Initial shape
        orbit = Square(side_length=4, color=BLUE, fill_opacity=0,stroke_width=polygon_stroke_width)
        target_dot = Dot(color=RED, radius=0.08).move_to(
            orbit.point_from_proportion(target_position.get_value())
        )

        #Rotation updater
        def orbit_updater_wrapper(mob):
            reference = mob.copy()
            def orbit_updater(inner_mob, dt):
                inner_mob.become(
                    reference.copy().rotate(
                        angle=orbit_angle.get_value() * 2 * np.pi * orbit_rotations_per_second,
                        about_point=ORIGIN # Changing this creates weird behavior
                    )
                )
            return orbit_updater

        orbit.add_updater(orbit_updater_wrapper(orbit))

        #Target dot updater
        def target_updater(mob, dt):
            mob.move_to(orbit.point_from_proportion(target_position.get_value()))

        target_dot.add_updater(target_updater)

        #Add initial objects to scene
        self.add(orbit, target_dot)

        #Traced line of the target dot
        target_tracer = TracedPath(lambda : target_dot.get_center(),stroke_width=1.5,stroke_color=WHITE,dissipating_time=tail_dissipation_time)

        dot_list = VGroup()
        self.add(target_tracer)

        #Iterative generations
        for gen in range(n_generations):

            # Get current vertices of the shape (evenly spaced for first generation)
            if gen == 0:
                current_vertices = [
                    target_dot.copy().move_to(orbit.point_from_proportion(min(p,1))).get_center()
                    for p in np.linspace(0, 1, n_vertices + 1)
                ][:-1]  # exclude last to avoid duplicating start. 
            else:
                current_vertices = new_coords  # from previous generation

            # Compute segment lengths and cumulative proportions
            segment_lengths = []
            total_length = 0
            n = len(current_vertices)
            for i in range(n):
                seg_len = np.linalg.norm(current_vertices[(i + 1) % n] - current_vertices[i])
                segment_lengths.append(seg_len)
                total_length += seg_len

            edge_proportions = [l / total_length for l in segment_lengths]
            cumulative_proportions = np.cumsum(edge_proportions)

            # Determine next n_vertices along perimeter starting AFTER 0
            # This ensures no duplicate of the starting vertex
            new_vertex_proportions = cumulative_proportions# This was "cumulative_proportions[:n_vertices]" earlier but this works aswell, because the linspace in the for gen in (...) loop filters the vertex proportions already
            print("New vertex proportions on iteration {} are: {}".format(gen+1,new_vertex_proportions))

            # Generate new coordinates based on the orbit and rotation. Also what the fuck is this list constructor abomination
            new_coords = [
                target_dot.copy().move_to(
                    orbit.point_from_proportion(min(1,p))
                ).rotate(
                    angle=p * 2 * np.pi * orbit_rotations_per_second,
                    about_point=ORIGIN
                ).get_center()
                for p in new_vertex_proportions
            ]

            # Optional: add reference dots for new vertices
            new_dots = VGroup([Dot(radius=new_dot_radius, color=WHITE).move_to(coord) for coord in new_coords])
            self.add(new_dots)
            
            self.play(
                Indicate(new_dots,scale_factor=1.4),
                FadeOut(dot_list[-1:]),
                run_time=new_dot_indication_time
            )

            # Remove old dots after FadeOut so the animation doesnt slow down while rendering a high amount of generations
            self.remove(dot_list[-1:])

            # Add the new dots to the list so they can be faded out in the next generation 
            dot_list.add(new_dots)

            # Animate rotation + target dot along perimeter
            self.play(
                orbit_angle.animate.set_value(1),
                target_position.animate.set_value(1),
                rate_func=linear,
                run_time=generation_runtime
            )

            # Morph current orbit to new polygon
            new_polygon = Polygon(*new_coords, color=BLUE, fill_opacity=0,stroke_width=polygon_stroke_width)
            self.play(
                Transform(orbit, new_polygon),
                Indicate(new_dots,scale_factor=1.4),
                run_time=transform_runtime
            )

            # Reset updater for new polygon to continue rotating. I feel like there should be a better way to do this but I haven't found it yet
            orbit.remove_updater(orbit.get_updaters())
            orbit.become(new_polygon)
            orbit.add_updater(orbit_updater_wrapper(orbit))

            # Reset trackers for next generation
            orbit_angle.set_value(0)
            target_position.set_value(0)



class MandelbrotZoom(Scene):

    def mandelbrot_image(self, width, height, x_center, y_center, zoom, max_iter=120, fname="frame.png"):
        """
        Computes the Mandelbrot set and saves it as a PNG that Manim can load.
        """
        scale = 1 / zoom
        x_min = x_center - 3.5 * scale
        x_max = x_center + 3.5 * scale
        y_min = y_center - 2.0 * scale
        y_max = y_center + 2.0 * scale

        xs = np.linspace(x_min, x_max, width)
        ys = np.linspace(y_min, y_max, height)

        mandelbrot = np.zeros((height, width))

        for iy, y in enumerate(ys):
            for ix, x in enumerate(xs):
                c = complex(x, y)
                z = 0
                n = 0
                while abs(z) <= 2 and n < max_iter:
                    z = z * z + c
                    n += 1
                mandelbrot[iy, ix] = n

        # Convert escape times → color image
        plt.imsave(fname, mandelbrot, cmap="magma")

        return fname

    def construct(self):

        resolution = 400
        max_iter = 150

        x_center = -0.75
        y_center = 0.0
        zoom = 1.0

        # Initial frame
        fname = self.mandelbrot_image(resolution, resolution, x_center, y_center, zoom, max_iter)
        img = ImageMobject(fname).set_height(6)
        self.add(img)

        # ZOOM ANIMATION
        def update_image(mob, alpha):
            z = 1 + alpha * 2000    # deep zoom 
            fname = f"frame_{int(alpha*1000)}.png"
            self.mandelbrot_image(resolution, resolution, x_center, y_center, z, max_iter, fname)
            new_img = ImageMobject(fname).set_height(6)
            mob.become(new_img)

        self.play(UpdateFromAlphaFunc(img, update_image), run_time=12, rate_func=smooth)
        self.wait()

        # optional clean up
        for f in os.listdir("."):
            if f.startswith("frame_") and f.endswith(".png"):
                os.remove(f)

class recursive_squares(MovingCameraScene):
    def construct(self):
        self.camera.frame_height = 15
        self.camera.frame_width = 15

        # Parameters
        initial_square_side_length = 4
        per_generation_rotations = 0.125 # 45 degrees in 2 pi rad
        per_generation_scale = np.sqrt(2)

        # Value Trackers
        generation_tracker = ValueTracker(0)
        generation_scale = ValueTracker(1)
        frame_rotation = ValueTracker(0)

        # Shape definitions
        i_square = Square(side_length=initial_square_side_length,color=GREEN,fill_opacity=0)
        #i_square_reference = i_square.copy().set_color(YELLOW)

        #self.add(i_square_reference)

        dot_1 = Dot(radius=0.07,color=WHITE)
        dot_2 = dot_1.copy()
        dot_3 = dot_1.copy()
        dot_4 = dot_1.copy()

        dot_group = VGroup(dot_1,dot_2,dot_3,dot_4)

        # Updaters
        def dots_updater(mob,dt):
            mob[0].move_to(i_square.point_from_proportion(generation_tracker.get_value()))
            mob[1].move_to(i_square.point_from_proportion(generation_tracker.get_value() + 0.25))
            mob[2].move_to(i_square.point_from_proportion(generation_tracker.get_value() + 0.5))
            mob[3].move_to(i_square.point_from_proportion(generation_tracker.get_value() + 0.75))

        dot_group.add_updater(dots_updater)

        def frame_updater_wrapper(object):
            reference = object.copy()
            def frame_updater(mob,dt):
                mob.become(reference.copy().rotate(angle=frame_rotation.get_value()*2*np.pi).scale(generation_scale.get_value()))

            return frame_updater

        i_square.add_updater(frame_updater_wrapper(i_square))

        # Always redraw lines
        def d1_to_d2():
            return Line(dot_1.get_center(),dot_2.get_center(),color=PINK)
        
        def d2_to_d3():
            return Line(dot_2.get_center(),dot_3.get_center(),color=PINK)
        
        def d3_to_d4():
            return Line(dot_3.get_center(),dot_4.get_center(),color=PINK)
        
        def d4_to_d1():
            return Line(dot_4.get_center(),dot_1.get_center(),color=PINK)
        
        line_12 = always_redraw(d1_to_d2)
        line_23 = always_redraw(d2_to_d3)
        line_34 = always_redraw(d3_to_d4)
        line_41 = always_redraw(d4_to_d1)
        


        # Add shapes
        self.add(i_square,dot_group)
        self.add(line_12,line_23,line_34,line_41)

        initial_scale = generation_scale.get_value()

        print("Pre generation scale = {}".format(initial_scale))

        self.play(
            generation_tracker.animate.set_value(per_generation_rotations),
            frame_rotation.animate.set_value(per_generation_rotations),
            generation_scale.animate.set_value(initial_scale+per_generation_scale),
            run_time=4
        )

        print("Post generation scale = {}".format(generation_scale.get_value()))

        print("ValueTracker values are generation_tracker = {}, frame_rotation = {}".format(generation_tracker.get_value(),frame_rotation.get_value()))

        # Add new dots for new inset square
        dot_11 = Dot(radius=0.07,color=RED)
        dot_21 = dot_11.copy()
        dot_31 = dot_11.copy()
        dot_41 = dot_11.copy()

        dot_group_2 = VGroup(dot_11,dot_21,dot_31,dot_41)

        print("Coordinates of new polygon vertices are = {}".format(
            [mobject.rotate(angle=generation_tracker.get_value()*2*np.pi).get_center() for mobject in dot_group]
        ))

        # Assign new polygon based on coordinates of where the previously moving dots were
        new_square = Polygon(
            *[mobject.rotate(angle=generation_tracker.get_value()*2*np.pi).get_center() for mobject in dot_group],
            #*[mobject.get_center() for mobject in dot_group],
            color=BLUE,
            fill_opacity = 0
        )

        def frame_updater_wrapper_1(object):
            reference = object.copy()
            def frame_updater_1(mob,dt):
                mob.become(reference.copy().rotate(angle=(frame_rotation.get_value()-per_generation_rotations)*2*np.pi).scale(generation_scale.get_value()))

            return frame_updater_1

        new_square.add_updater(frame_updater_wrapper_1(new_square))


        def dots_updater_1_wrapper(object_of_proportion):
            def dots_updater_1(mob,dt):
                mob[0].move_to(object_of_proportion.point_from_proportion(generation_tracker.get_value()))
                mob[1].move_to(object_of_proportion.point_from_proportion(generation_tracker.get_value() + 0.25))
                mob[2].move_to(object_of_proportion.point_from_proportion(generation_tracker.get_value() + 0.5))
                mob[3].move_to(object_of_proportion.point_from_proportion(generation_tracker.get_value() + 0.75))

            return dots_updater_1

        dot_group_2.add_updater(dots_updater_1_wrapper(new_square))

        # New always redraw lines
        def d1_to_d2_1():
            return Line(dot_11.get_center(),dot_21.get_center(),color=WHITE)
        
        def d2_to_d3_1():
            return Line(dot_21.get_center(),dot_31.get_center(),color=WHITE)
        
        def d3_to_d4_1():
            return Line(dot_31.get_center(),dot_41.get_center(),color=WHITE)
        
        def d4_to_d1_1():
            return Line(dot_41.get_center(),dot_11.get_center(),color=WHITE)
        

        line_12_1 = always_redraw(d1_to_d2_1)
        line_23_1 = always_redraw(d2_to_d3_1)
        line_34_1 = always_redraw(d3_to_d4_1)
        line_41_1 = always_redraw(d4_to_d1_1)

        self.add(new_square,dot_group_2,line_12_1,line_23_1,line_34_1,line_41_1)

        

        # Play rotation transform again
        self.play(
            generation_tracker.animate.set_value(2/8),
            frame_rotation.animate.set_value(2/8),
            run_time=4
        )

        self.wait(10)

class nostalrium(MovingCameraScene):
    def construct(self):
        self.camera.frame_height = 20
        self.camera.frame_width = 20
        
        self.t_offset = 0
        
        def time_tracker(dt):
            self.t_offset += dt
        self.add_updater(time_tracker)


        # Parameters
        group1_radial_offset = 2
        group2_radial_offset = 4.5

        # Shape definitions

        # Define radial offset
        def shift_radially(mob, dr, about_point=ORIGIN):
            # Current vector from about_point to mob
            v = mob.get_center_of_mass() - about_point
            r, theta = np.linalg.norm(v), np.arctan2(v[1], v[0])
            # New position with increased radius
            new_pos = about_point + (r + dr) * np.array([np.cos(theta), np.sin(theta), 0])
            mob.move_to(new_pos)
            return mob
        
        center_square = Square(side_length=1,color=WHITE,fill_opacity=0)
        tri_r = Polygon(
            [0.5,0.5,0],
            [1.5,0,0],
            [0.5,-0.5,0],
            color=WHITE,
            fill_opacity=0
        )
        

        tri_t = tri_r.copy().rotate(angle=90*DEGREES,about_point=ORIGIN)
        tri_l = tri_r.copy().rotate(angle=180*DEGREES,about_point=ORIGIN)
        tri_b = tri_r.copy().rotate(angle=270*DEGREES,about_point=ORIGIN)

        center_group_angle_tracker = ValueTracker(0)
        center_group_radial_offset_tracker = ValueTracker(0.5)

        outside_group_1_angle_tracker = ValueTracker(0)
        outside_group_1_radial_offset_tracker = ValueTracker(0)


        # Updaters for members of center group
        def tri_r_updater_wrapper(object):
            reference = object.copy()
            def center_triangle_updater(mob,dt):
                mob.become(shift_radially(reference.copy().rotate(angle=center_group_angle_tracker.get_value(),about_point=ORIGIN),center_group_radial_offset_tracker.get_value()))
            return center_triangle_updater

        def tri_t_updater_wrapper(object):
            reference = object.copy()
            def center_triangle_updater(mob,dt):
                mob.become(shift_radially(reference.copy().rotate(angle=center_group_angle_tracker.get_value()+np.pi/2,about_point=ORIGIN),center_group_radial_offset_tracker.get_value()))
            return center_triangle_updater

        def tri_l_updater_wrapper(object):
            reference = object.copy()
            def center_triangle_updater(mob,dt):
                mob.become(shift_radially(reference.copy().rotate(angle=center_group_angle_tracker.get_value()+np.pi,about_point=ORIGIN),center_group_radial_offset_tracker.get_value()))
            return center_triangle_updater

        def tri_b_updater_wrapper(object):
            reference = object.copy()
            def center_triangle_updater(mob,dt):
                mob.become(shift_radially(reference.copy().rotate(angle=center_group_angle_tracker.get_value()+np.pi*1.5,about_point=ORIGIN),center_group_radial_offset_tracker.get_value()))
            return center_triangle_updater

        tri_r.add_updater(tri_r_updater_wrapper(tri_r))
        tri_t.add_updater(tri_r_updater_wrapper(tri_t))
        tri_l.add_updater(tri_r_updater_wrapper(tri_l))
        tri_b.add_updater(tri_r_updater_wrapper(tri_b))

        center_group = VGroup(center_square,tri_r,tri_t,tri_l,tri_b)

        # Updaters for members of outer group1
        def tri_our_updater_wrapper(object):
            reference = object.copy()
            def outer_triangle_updater(mob,dt):
                mob.become(shift_radially(reference.copy().rotate(angle=outside_group_1_angle_tracker.get_value(),about_point=ORIGIN),outside_group_1_radial_offset_tracker.get_value()))
            return outer_triangle_updater

        def tri_oul_updater_wrapper(object):
            reference = object.copy()
            def outer_triangle_updater(mob,dt):
                mob.become(shift_radially(reference.copy().rotate(angle=outside_group_1_angle_tracker.get_value(),about_point=ORIGIN),outside_group_1_radial_offset_tracker.get_value()))
            return outer_triangle_updater

        def tri_odl_updater_wrapper(object):
            reference = object.copy()
            def outer_triangle_updater(mob,dt):
                mob.become(shift_radially(reference.copy().rotate(angle=outside_group_1_angle_tracker.get_value(),about_point=ORIGIN),outside_group_1_radial_offset_tracker.get_value()))
            return outer_triangle_updater

        def tri_odr_updater_wrapper(object):
            reference = object.copy()
            def outer_triangle_updater(mob,dt):
                mob.become(shift_radially(reference.copy().rotate(angle=outside_group_1_angle_tracker.get_value(),about_point=ORIGIN),outside_group_1_radial_offset_tracker.get_value()))
            return outer_triangle_updater
        


        base_center_triangle = Polygon(
            [-0.5,0,0],
            [0.5,0,0],
            [0,1,0],
            color=WHITE,
            fill_opacity=0
        ).shift(DOWN*0.5)

        tri_our = base_center_triangle.copy().rotate(angle=135*DEGREES).shift(UP*group1_radial_offset+RIGHT*group1_radial_offset)
        tri_oul = base_center_triangle.copy().rotate(angle=215*DEGREES).shift(UP*group1_radial_offset+LEFT*group1_radial_offset)
        tri_odl = base_center_triangle.copy().rotate(angle=305*DEGREES).shift(DOWN*group1_radial_offset+LEFT*group1_radial_offset)
        tri_odr = base_center_triangle.copy().rotate(angle=45*DEGREES).shift(DOWN*group1_radial_offset+RIGHT*group1_radial_offset)

        tri_our.add_updater(tri_our_updater_wrapper(tri_our))
        tri_oul.add_updater(tri_oul_updater_wrapper(tri_oul))
        tri_odl.add_updater(tri_odl_updater_wrapper(tri_odl))
        tri_odr.add_updater(tri_odr_updater_wrapper(tri_odr))

        outside_triangle_group_1 = VGroup(tri_our,tri_oul,tri_odl,tri_odr)

        tri_oor = base_center_triangle.copy().rotate(angle=90*DEGREES).shift(RIGHT*group2_radial_offset)
        tri_oot = base_center_triangle.copy().rotate(angle=180*DEGREES).shift(UP*group2_radial_offset)
        tri_ool = base_center_triangle.copy().rotate(angle=270*DEGREES).shift(LEFT*group2_radial_offset)
        tri_ood = base_center_triangle.copy().rotate(angle=360*DEGREES).shift(DOWN*group2_radial_offset)

        outside_triangle_group_2 = VGroup(tri_oor,tri_oot,tri_ool,tri_ood)

        self.add(center_group,outside_triangle_group_1,outside_triangle_group_2)
        self.play(
            #center_group_angle_tracker.animate.set_value(2*np.pi),
            Rotate(center_group,angle=360*DEGREES),
            Rotate(outside_triangle_group_1,angle=360*DEGREES),
            center_group_radial_offset_tracker.animate.set_value(2),
            outside_group_1_radial_offset_tracker.animate.set_value(-5),
            run_time=4,
            rate_func = linear
        )

class sympy_square(MovingCameraScene):
    def construct(self):

        # Parameters
        ray_starting_position = [-0.5,0.5,0]
        ray_starting_direction_nonnormalised = [1.5,1,0]
        ray_head_speed = 5 # Units per second


        # Preformatting some stuff like normalising vectors for directions
        ray_starting_direction = np.linalg.norm(ray_starting_direction_nonnormalised)**(-1) * np.array(ray_starting_direction_nonnormalised)

        # Shape definitions

        ray_head = Dot(ray_starting_position,radius=0.05,color=RED)
        ray_tracer = TracedPath(lambda: ray_head.get_center(),color=WHITE,stroke_width=1,dissipating_time=8)
        self.add(ray_head,ray_tracer)

        bounding_box_coords = np.array([
            [-2,-2,0],
            [-2,1,0],
            [2,3,0],
            [2,-5,0]
        ])

        #bounding_box_coords = np.array([
        #    [ 1.22464680e-16, -2.00000000e+00,  0.00000000e+00],
        #    [ 5.87785252e-01, -8.09016994e-01,  0.00000000e+00],
        #    [ 1.90211303e+00, -6.18033989e-01,  0.00000000e+00],
        #    [ 9.51056516e-01,  3.09016994e-01,  0.00000000e+00],
        #    [ 1.17557050e+00,  1.61803399e+00,  0.00000000e+00],
        #    [ 6.12323400e-17,  1.00000000e+00,  0.00000000e+00],
        #    [-1.17557050e+00,  1.61803399e+00,  0.00000000e+00],
        #    [-9.51056516e-01,  3.09016994e-01,  0.00000000e+00],
        #    [-1.90211303e+00, -6.18033989e-01,  0.00000000e+00],
        #    [-5.87785252e-01, -8.09016994e-01,  0.00000000e+00]]
        #)

        bounding_box = Polygon(
            *bounding_box_coords,
            color=WHITE,
            fill_opacity=0
        )
        self.add(bounding_box)

        # Array of vectors connecting the vertices in order of Polygon creation
        vertex_connectors = np.zeros_like(bounding_box_coords)

        for i in range(len(bounding_box_coords)):
            vertex_connectors[i] = bounding_box_coords[(i+1) % len(bounding_box_coords)]-bounding_box_coords[i]

        def line_intersection(ray_head_coordinate,ray_direction,vertex_coordinate,vertex_connector_direction,eps=1e-8):
            t, s = sp.symbols(' t s', real=True)

            equation = sp.Matrix(ray_head_coordinate) + t * sp.Matrix(ray_direction) - sp.Matrix(vertex_coordinate + s*vertex_connector_direction)

            solution = sp.solve([equation[0], equation[1]], [t, s], dict=True)

            if not solution:
                return None # Either somethings parallel or there is no intersection (can only happen in 3D but wont apply here, parallel is biggest issue)
            
            t_value = float(solution[0][t])
            s_value = float(solution[0][s])

            if (t_value > eps) and (s_value >= -eps) and (s_value <= 1+ eps):
                return ray_head_coordinate + t_value * ray_direction
            else:
                return ray_head_coordinate + t_value * ray_direction * 1e3
        
        # Compute all vertex connector intersections
        self.intersections = []

        def get_line_intersections(vertex_coordinates,vertex_connectors,ray_head_coordinate,ray_direction):
            self.intersections = []
            for i in range(len(vertex_coordinates)):
                current_line_intersection = line_intersection(ray_head_coordinate,ray_direction,vertex_coordinates[i],vertex_connectors[i])
                self.intersections.append(current_line_intersection) 

        self.intersection_distances = []
        self.new_ray_head_line = []
        self.shortest_intersection_index = []

        def find_closest_intersection(ray_head_coordinate,intersection_coordinates):
            self.intersection_distances = []
            self.new_ray_head_line = []
            self.shortest_intersection_index = []
            for i in range(len(intersection_coordinates)):
                print("Ray head coordinate when find_closest_intersection is called {}".format(ray_head_coordinate))
                print("Current iteration {} of intersection coordinate: {}".format(i,intersection_coordinates[i]))
                current_distance = np.linalg.norm(ray_head_coordinate-intersection_coordinates[i])
                if current_distance > 0.01:
                    self.intersection_distances.append(current_distance)
                else:
                    # Add a huge artificial distance so the point doesnt reflect with "itself" to a very close intersection between new direction and the edge it just hit
                    self.intersection_distances.append(1e3)
            
            self.shortest_intersection_index = np.argmin(self.intersection_distances)
            print("Shortest intersection is at: {}".format(self.intersections[self.shortest_intersection_index]))
            print("Distance to closest intersection is {}".format(self.intersection_distances[self.shortest_intersection_index]))
            self.new_ray_head_line = self.intersections[self.shortest_intersection_index] - ray_head_coordinate

        def move_to_intersection(mob):
            self.play(
                mob.animate.move_to(mob.get_center()+self.new_ray_head_line),
                run_time= 1/ray_head_speed * np.linalg.norm(self.new_ray_head_line),
                rate_func=linear
            )

        self.new_velocity_direction = []

        def get_new_velocity_direction(ray_direction):
            self.new_velocity_direction = []
            polygon_tangent = vertex_connectors[self.shortest_intersection_index] / np.linalg.norm(vertex_connectors[self.shortest_intersection_index])
            polygon_tangent_normal = np.array([-polygon_tangent[1],polygon_tangent[0],0])
            scalar_product = np.dot(ray_direction,polygon_tangent_normal)
            self.new_velocity_direction = ray_direction-2*scalar_product*polygon_tangent_normal


        get_line_intersections(bounding_box_coords, vertex_connectors,ray_head.get_center(),ray_starting_direction)
        find_closest_intersection(ray_head.get_center(),self.intersections)
        move_to_intersection(ray_head)
        get_new_velocity_direction(ray_starting_direction)
        print("New ray head direction is: {}".format(self.new_ray_head_line))

        #self.play(
        #    ray_head.animate.move_to(ray_head.get_center()+self.new_ray_head_direction),
        #    run_time=2
        #)

        reflections_to_render = 30

        for i in range(reflections_to_render):
            get_line_intersections(bounding_box_coords, vertex_connectors,ray_head.get_center(),self.new_velocity_direction)
            find_closest_intersection(ray_head.get_center(),self.intersections)
            move_to_intersection(ray_head)
            get_new_velocity_direction(self.new_ray_head_line/np.linalg.norm(self.new_ray_head_line))

class sympy_square_parameterpath(MovingCameraScene):
    def construct(self):
        self.camera.frame_height = 20
        self.camera.frame_width = 20

        # Parameters
        ray_starting_position = [0,0,0]
        ray_starting_direction_nonnormalised = [1,0.25,0]
        ray_head_speed = 5 # Units per second
        ray_tail_dissipation_time = 0.2
        ray_tracer_stroke_width = 1.2
        bounding_box_stroke_width = 1.95
        animation_run_time = 100
        reflections_to_render = 103

        path_position = ValueTracker(0)


        # Preformatting some stuff like normalising vectors for directions
        ray_starting_direction = np.linalg.norm(ray_starting_direction_nonnormalised)**(-1) * np.array(ray_starting_direction_nonnormalised)

        # Shape definitions

        ray_head = Dot(ray_starting_position,radius=0.035,color=RED)
        ray_tracer = TracedPath(lambda: ray_head.get_center(),color=WHITE,stroke_width=ray_tracer_stroke_width,dissipating_time=ray_tail_dissipation_time)
        self.add(ray_head,ray_tracer)


        def load_to_3d(coords):
            return np.column_stack([coords, np.zeros(len(coords))])
        bounding_box_coords = load_to_3d(np.loadtxt("26112025_141209.txt")) 
        #print(type(bounding_box_coords))


        #bounding_box_coords = np.array([
        #    [-0.5,-0.5,0],
        #    [0,0.5,0],
        #    [0.5,-0.5,0]
        #])

        #bounding_box_coords = np.array([
        #    [-2,-2,0],
        #    [-2,1,0],
        #    [2,3,0],
        #    [2,-5,0]
        #])

        #bounding_box_coords = np.array([
        #    [ 1.22464680e-16, -2.00000000e+00,  0.00000000e+00],
        #    [ 5.87785252e-01, -8.09016994e-01,  0.00000000e+00],
        #    [ 1.90211303e+00, -6.18033989e-01,  0.00000000e+00],
        #    [ 9.51056516e-01,  3.09016994e-01,  0.00000000e+00],
        #    [ 1.17557050e+00,  1.61803399e+00,  0.00000000e+00],
        #    [ 6.12323400e-17,  1.00000000e+00,  0.00000000e+00],
        #    [-1.17557050e+00,  1.61803399e+00,  0.00000000e+00],
        #    [-9.51056516e-01,  3.09016994e-01,  0.00000000e+00],
        #    [-1.90211303e+00, -6.18033989e-01,  0.00000000e+00],
        #    [-5.87785252e-01, -8.09016994e-01,  0.00000000e+00]]
        #)

        bounding_box = Polygon(
            *bounding_box_coords,
            color=WHITE,
            stroke_width=bounding_box_stroke_width,
            fill_opacity=0
        )
        self.add(bounding_box)

        # Array of vectors connecting the vertices in order of Polygon creation
        vertex_connectors = np.zeros_like(bounding_box_coords)

        for i in range(len(bounding_box_coords)):
            vertex_connectors[i] = bounding_box_coords[(i+1) % len(bounding_box_coords)]-bounding_box_coords[i]

        def line_intersection(ray_head_coordinate,ray_direction,vertex_coordinate,vertex_connector_direction,eps=1e-8):
            t, s = sp.symbols(' t s', real=True)

            equation = sp.Matrix(ray_head_coordinate) + t * sp.Matrix(ray_direction) - sp.Matrix(vertex_coordinate + s*vertex_connector_direction)

            solution = sp.solve([equation[0], equation[1]], [t, s], dict=True)

            if not solution:
                #return None # Either somethings parallel or there is no intersection (can only happen in 3D but wont apply here, parallel is biggest issue)
                return ray_head_coordinate + 1e6 * ray_direction
            
            t_value = float(solution[0][t])
            s_value = float(solution[0][s])

            if (t_value > eps) and (s_value >= -eps) and (s_value <= 1+ eps):
                return ray_head_coordinate + t_value * ray_direction
            else:
                return ray_head_coordinate + t_value * ray_direction * 1e3
        
        # Compute all vertex connector intersections
        self.intersections = []

        def get_line_intersections(vertex_coordinates,vertex_connectors,ray_head_coordinate,ray_direction):
            self.intersections = []
            for i in range(len(vertex_coordinates)):
                current_line_intersection = line_intersection(ray_head_coordinate,ray_direction,vertex_coordinates[i],vertex_connectors[i])
                self.intersections.append(current_line_intersection) 

        self.intersection_distances = []
        self.new_ray_head_line = []
        self.shortest_intersection_index = []

        def find_closest_intersection(ray_head_coordinate,intersection_coordinates):
            self.intersection_distances = []
            self.new_ray_head_line = []
            self.shortest_intersection_index = []
            for i in range(len(intersection_coordinates)):
                #print("Ray head coordinate when find_closest_intersection is called {}".format(ray_head_coordinate))
                #print("Current iteration {} of intersection coordinate: {}".format(i,intersection_coordinates[i]))
                current_distance = np.linalg.norm(ray_head_coordinate-intersection_coordinates[i])
                if current_distance > 0.01:
                    self.intersection_distances.append(current_distance)
                else:
                    # Add a huge artificial distance so the point doesnt reflect with "itself" to a very close intersection between new direction and the edge it just hit
                    self.intersection_distances.append(1e3)
            
            self.shortest_intersection_index = np.argmin(self.intersection_distances)
            #print("Shortest intersection is at: {}".format(self.intersections[self.shortest_intersection_index]))
            print("Distance to closest intersection is {}".format(self.intersection_distances[self.shortest_intersection_index]))
            self.new_ray_head_line = self.intersections[self.shortest_intersection_index] - ray_head_coordinate

        def move_to_intersection(mob):
            self.play(
                mob.animate.move_to(mob.get_center()+self.new_ray_head_line),
                run_time= 1/ray_head_speed * np.linalg.norm(self.new_ray_head_line),
                rate_func=linear
            )

        self.new_velocity_direction = []

        def get_new_velocity_direction(ray_direction):
            self.new_velocity_direction = []
            polygon_tangent = vertex_connectors[self.shortest_intersection_index] / np.linalg.norm(vertex_connectors[self.shortest_intersection_index])
            polygon_tangent_normal = np.array([-polygon_tangent[1],polygon_tangent[0],0])
            scalar_product = np.dot(ray_direction,polygon_tangent_normal)
            self.new_velocity_direction = ray_direction-2*scalar_product*polygon_tangent_normal

        ray_path_coordinates = [ray_starting_position.copy()]
        
        get_line_intersections(bounding_box_coords, vertex_connectors,ray_head.get_center(),ray_starting_direction)
        find_closest_intersection(ray_head.get_center(),self.intersections)
        ray_path_coordinates.append(ray_path_coordinates[-1]+self.new_ray_head_line)
        #move_to_intersection(ray_head)
        get_new_velocity_direction(ray_starting_direction)
        #print("New ray head direction is: {}".format(self.new_ray_head_line))


        
        #ray_path_coordinates = np.zeros([reflections_to_render+1,3])
        #ray_path_coordinates = [ray_starting_position.copy()]
        #ray_path_coordinates[0,:] = ray_starting_position

        for i in range(reflections_to_render):
            print("Iteration {} / {}".format(i,reflections_to_render))
            get_line_intersections(bounding_box_coords, vertex_connectors,ray_path_coordinates[-1],self.new_velocity_direction)
            find_closest_intersection(ray_path_coordinates[-1],self.intersections)
            new_point = ray_path_coordinates[-1] + self.new_ray_head_line
            ray_path_coordinates.append(new_point)
            #move_to_intersection(ray_head)
            get_new_velocity_direction(self.new_ray_head_line/np.linalg.norm(self.new_ray_head_line))

        # Construct Polygon along path of ray
        ray_path = VMobject()
        ray_path.set_points_as_corners(ray_path_coordinates)

        # Updater for ray_head
        def ray_head_updater(mob,dt):
            mob.move_to(ray_path.point_from_proportion(path_position.get_value()))

        ray_head.add_updater(ray_head_updater)
        self.play(
            path_position.animate.set_value(1),
            run_time=animation_run_time,
            rate_func=linear
        )


