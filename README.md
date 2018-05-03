Function Fitting3DPointsToLine
To fitting 3d points, which are all on one plane, into a line



Input: txt file (testfile*.txt)
       first row is the normal of the plane on which all points are
       every row of the rest are coordinates (xyz) of these points
       
    
    
      
Output: two points(head & tail) on the fitting line



flag_head_tail == 0: head and tail are the projection of the first two points in points_list
flag_head_tail == 1: head and tail are the projection of the head and tail points in points_list



fitting_mode == 0: 2d PCA
fitting_mode == 1: Linear Regression(Perpendicular least squares fitting)
fitting_mode == 2: Linear Regression(Vertical least squares fitting)
