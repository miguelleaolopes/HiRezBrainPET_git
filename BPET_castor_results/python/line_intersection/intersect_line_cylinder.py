import numpy as np
from numpy import sqrt
import matplotlib.pyplot as plt


def intersect_line_cylinder(point, direction, radius, height, show=False, verbose=False):
    """
    Find the intersection points of a line and a cylinder.

    Parameters
    ----------
    point : list or np.array (x, y, z)
        Coordinates of a point on the line.
    direction : list or np.array (vx, vy, vz)
        Direction vector of the line. (if not a unit vector, it will be normalized)
    radius : float
        Radius of the cylinder.
    height : float
        Height of the cylinder.
    show : bool, optional
        Flag to show the line and cylinder. Default is False.
    verbose : bool, optional
        Flag to print messages. Default is False.


    Returns
    -------
    array or None
        If the line intersects the cylinder, returns an array of two points, each represented as a list of three coordinates.
        If the line does not intersect the cylinder, returns None.
    """
    # Ensure direction is a unit vector
    direction = direction/np.linalg.norm(direction)
    
    # Vector from point to cylinder center
    cylinder_center = np.array([0, 0, -height/2])
    cylinder_center_to_point = point - cylinder_center
    z_direction = np.array([0, 0, 1])

    # Calculate the coefficients of the quadratic equation
    a = np.dot(direction, direction)-np.dot(direction, z_direction)**2
    b = 2*(np.dot(cylinder_center_to_point, direction) - np.dot(cylinder_center_to_point, z_direction)*np.dot(direction, z_direction))
    c = np.dot(cylinder_center_to_point, cylinder_center_to_point) - np.dot(cylinder_center_to_point, z_direction)**2 - radius**2

    # Calculate the discriminant
    discriminant = b**2 - 4*a*c

    # If discriminant is negative, there are no real roots
    if discriminant < 0 :
        f2Points = False
        if verbose:
            print('Line does not intersect cylinder')    # If discriminant is zero, there is one real root
    elif discriminant == 0:
        f2Points = False
        # Calculate the intersection point
        t = -b/(2*a)
        intersection_point = point + t*direction
        if verbose:
            print('Line intersects cylinder at only 1 point:', intersection_point)
        # If discriminant is positive, there are two real roots
    elif np.isnan(discriminant):
        f2Points = False
        if verbose:
            print('Line does not intersect cylinder because discriminant is NaN')
    else:
        f2Points = True
        # Calculate the intersection points
        t1 = (-b + sqrt(discriminant))/(2*a)
        t2 = (-b - sqrt(discriminant))/(2*a)
        intersection_point1 = point + t1*direction
        intersection_point2 = point + t2*direction

    ######################################################################
    # Check if intersection points are within the height of the cylinder
    if f2Points:
        if intersection_point1[2] < -height/2 or intersection_point1[2] > height/2 or intersection_point2[2] < -height/2 or intersection_point2[2] > height/2:
            if verbose:
                print('Line is not within the height of the cylinder')
            f2Points = False
    ######################################################################

    # Show line and cylinder
    if show:
        # Number of points to draw the line in 3d
        npts = 1000

        # Plot Cylinder (surface)
        theta = np.linspace(0, 2*np.pi, npts)
        x = radius*np.cos(theta)
        y = radius*np.sin(theta)
        z = np.full_like(theta, height/2)
        ax = plt.axes(projection='3d')
        ax.plot3D(x, y, z, 'b', label='Cylinder')
        ax.plot3D(x, y, -z, 'b')
        # Plot Cylinder (edges)
        ax.plot3D([radius, radius], [0, 0], [height/2, -height/2], 'b')
        ax.plot3D([-radius, -radius], [0, 0], [height/2, -height/2], 'b')
        ax.plot3D([0, 0], [radius, radius], [height/2, -height/2], 'b')
        ax.plot3D([0, 0], [-radius, -radius], [height/2, -height/2], 'b')

        # Plot line (3D)
        t = np.linspace(-radius*1.5, radius*1.5, npts)
        x = point[0] + t*direction[0]
        y = point[1] + t*direction[1]
        z = point[2] + t*direction[2]
        ax.plot3D(x, y, z, 'r', label='Line')

        # Plot intersection points
        if f2Points:
            ax.scatter3D(intersection_point1[0], intersection_point1[1], intersection_point1[2], color='purple', label='Intersection Points 1')
            ax.scatter3D(intersection_point2[0], intersection_point2[1], intersection_point2[2], color='g', label='Intersection Points 2')

        # Plot settings
        ax.legend(loc='best')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('Line-Cylinder Intersection')
        plt.show()


    # Return intersection points (if they exist)
    if f2Points: return np.array([intersection_point1, intersection_point2])
    else: return None