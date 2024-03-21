import numpy as np
from numpy import sqrt
import matplotlib.pyplot as plt


def intersect_line_circle(A,B,C,r,show=False, verbose=False):
    """
    Find the intersection points of a line and a circle.

    Parameters
    ----------
    A : float
        Coefficient of the line equation Ax + By + C = 0.
    B : float
        Coefficient of the line equation Ax + By + C = 0.
    C : float
        Coefficient of the line equation Ax + By + C = 0.
    r : float
        Radius of the circle.
    show : bool, optional
        Flag to show the line and circle. Default is False.
    verbose : bool, optional
        Flag to print messages. Default is False.


    Returns
    -------
    list or int
        If the line intersects the circle, returns a list of two points, each represented as a list of two coordinates.
        If the line does not intersect the circle, returns None.

    Notes
    -----
    This function uses the method described in https://cp-algorithms.com/geometry/circle-line-intersection.html.
    It assumes the circle is centered at (0, 0). If it is not, the line and circle are translated to (0, 0), the intersection points are calculated, and then the result is translated back to the original coordinates.

    """
    
    # Denominator for the quadratic equation
    denom = A*A+B*B 

    # Coordinates of closest point to circle center
    x0 = - A*C/denom
    y0 = - B*C/denom
    # Distance from (x0,y0) to circle center
    d0 = sqrt(x0*x0 + y0*y0)

    # Intersection points
    if d0 < r: # Check if line intersects circle (2 points)
        # Value to add/subtract for reaching intersection points
        d   = r*r - C*C/denom # d^2 = r^2 - C^2/denom
        m   = sqrt(d/denom)   # m^2 = d^2/denom

        # Coords of intersection points (ax,ay) and (bx,by)
        ax  = x0 + B*m 
        ay  = y0 - A*m
        bx  = x0 - B*m
        by  = y0 + A*m
        f2Points = True
    else:
        f2Points = False
        if verbose:
            print('Line does not intersect circle.')

    # Show line and circle
    if show:
        # Number of points to draw the line
        npts = 1000

        # Plot circle
        theta = np.linspace(0, 2*np.pi, npts)
        circle_x = r*np.cos(theta)
        circle_y = r*np.sin(theta)
        plt.plot(circle_x, circle_y, 'b-', label='Circle')

        # Plot line (with a little more space) with intersection points (if they exist)
        x = np.linspace(-r*1.1, r*1.1, npts)
        y = (-A*x - C)/B
        # Close to a vertical line
        if y.max() > 1.1*r or y.min() < -1.1*r:
            # Cut high values
            idx = np.where(y > 1.1*r)
            x[idx] = np.nan
            y[idx] = np.nan 
            # Cut low values
            idx = np.where(y < -1.1*r)
            x[idx] = np.nan
            y[idx] = np.nan
        
        plt.plot(x, y, 'r-', label='Line')
        if f2Points:
            plt.plot(ax, ay, 'go', label='Intersection points')
            plt.plot(bx, by, 'go')

        # Plot settings
        plt.axis('equal')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Line-Circle Intersection')
        plt.legend(loc='best')
        plt.grid()
        plt.show() 

    # Return intersection points (if they exist)
    if f2Points: return [[ax,ay],[bx,by]]
    else: return None