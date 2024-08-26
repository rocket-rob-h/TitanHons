-- Make  the  simulation axisymetric:

config.axisymmetric = true

-- Variables

NoseRadius = 0.02615   -- [m]              (radius of  the spherecone) 90 mm model
CapsuleRadius  = 0.09   -- [m]              (radius  (frontal)  of the capsule)
ConeTheta =  60     --  [degrees]  (spherecone angle)
Offset = 0.005 -- [m] (initial offset from geometry of left inf low boundary)

function linearPart(x,  m, c)
--  """Just  a linear  relation"""
return  m*x  + c
end

function inverseLinearPart(y, m, c)
--  """Just  an  inverse linear relation"""
return  (y -  c)/m
end



--  The  coefficients of  the linear  part
m =  math.tan(ConeTheta *  math.pi  / 180)  --  gradient of linear  part
c =  (math.sqrt(1+m^2)-m)*NoseRadius         --  y-intercept of linear  part


--Coordinates
Bx = ( (2*NoseRadius - 2*m*c) + math.sqrt( (2*m*c - 2*NoseRadius)^2 - 4*(1+m^2)*(c^2) ) ) / (2 * (1+m^2))
By =  linearPart(Bx, m,  c)
Cy = 0.2
Cx = inverseLinearPart(Cy, m, c)
A = Vector3:new{x = 0.0, y = 0.0}
B =  Vector3:new{x = Bx, y  = By}    --  the  intersection  of the line  and  arc
C =  Vector3:new{x  = Cx,  y   =Cy}
CircCentre  = Vector3:new{x = NoseRadius,  y = 0.0}
D = C  + Vector3:new{x = 0.0, y =  Offset}
E  = A  -  Vector3:new{x = Offset,  y = 0.0}
bez_p1  = E + Vector3:new{x  = 0.0,  y = 2*Offset}
bez_p2  = D -  Vector3:new{x  = 2*Offset,  y = Offset}

-- Define the lines, arc,  and northern  Bezier:
--Lines, arcs
AB  = Arc:new{p0 = A,  centre = CircCentre, p1  = B}
BC = Line:new{p0  = B, p1  = C}
DC =  Line:new{p0  = D, p1 = C}
ED =  Bezier:new{points  = {E,  bez_p1,  bez_p2,  D}}
EA  =  Line:new{p0  = E,  p1 = A}

-- Consolidate the lines AB and BC into AC.  This is done so that the east boundary of the fluid
--block  can  be given  as  one line.

--Lines, arcs
AB  = Arc:new{p0 = A,  centre = CircCentre, p1  = B}
BC = Line:new{p0  = B, p1  = C}
DC =  Line:new{p0  = D, p1 = C}
ED =  Bezier:new{points  = {E,  bez_p1,  bez_p2,  D}}
EA  =  Line:new{p0  = E,  p1 = A}
AC = Polyline:new{segments  = {AB,BC}}

-- Finally, the surface, grid, fluid block array can be built. Here an example of clustering (toward
--the shock and stagnation line) is demonstrated. A shock-fitting boundary condition is set on
--the west edge of the fluid block so that the shock can be captured when the simulation is run.

niv  =32+1; njv  =128+1;
cluster_generic  = RobertsFunction:new{end0=true, end1=false, beta  =  1.01}
clusterlist =  {north=none,  south=none, east=cluster_generic, west=cluster_generic}
grid  = StructuredGrid:new{psurface=surface,  niv=niv,  njv=njv,  cfList=clusterlist}
blocks  = FBArray:new{grid=grid,
initialState=initial,
bcList={
north=OutFlowBC_Simple:new{},
west=InFlowBC_ShockFitting:new{f lowState=inf low}},
nib=6, njb=8}