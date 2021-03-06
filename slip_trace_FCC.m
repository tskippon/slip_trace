function [ plane, direction, traceMismatch,slipDirAngle] = slip_trace_FCC(imgfile, ebsdfile,n)
%SLIP_TRACE Takes in an SEM image file and EBSD data from a sample
%containing slip traces and identifies the slip traces selected graphically
%by the user.  The n most likely slip systems are returned 
%(if n is unspecified it will return the 3 most likely).

%   The slip plane is identified by comparing the direction formed by the
%   slip trace to the slip traces tha would be produced by the slip planes
%   known to be active in Zr.

% Returns: Listed below are the variables that get returned when the code
% is run.

%         plane - the n most likely slip planes

%         direction - the slip directions

%         traceMismatch - the angle (in degrees) between the user selected
%         slip trace and theoretical slip trace based on crystal symmetry

%         slipDirAngle - the angle (in degrees) that the slip direction
%         points out of the SEM image plane.  If this angle is 0 then this
%         slip system cannot produce visible slip traces in the image.

% Example: Below is how to call the function.

% [plane,direction,traceMismatch,slipDirAngle]=slip_trace_FCC('imagefile.tif','ebsdmap.ctf',3)


%Set n to 3 if not entered by the user.
if(nargin<3)
    n=3;
end

%Read in image file and have user click on slip trace
image=imread(imgfile);
figure(1);
imshow(image);

uiwait(msgbox('Click two points on a single slip trace.','Select Corners','modal'));
c1 = ginput(1);
c2 = ginput(1);
close(1);
line=normalize(vector3d(c2(1)-c1(1),c2(2)-c1(2),0));



%Read in ebsd file and calculate grains
CS = {'notIndexed',crystalSymmetry('m-3m', [3.524 3.524 3.524], 'mineral', 'Nickel', 'color', 'light blue')};
ebsd = loadEBSD(ebsdfile,CS,'interface','ctf','convertEuler2SpatialReferenceFrame');
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',5*degree);

%Display plot of grains and have user select the one that their SEM image
%is located in.
figure(2)
plot(grains(grains.phase==1),grains(grains.phase==1).meanOrientation);
uiwait(msgbox('Click on the grain that is being indented to get crystal orientation.   Make sure SEM image and EBSD map are oriented the same way!','choose grain','modal'));
g=ginput(1);
close 2;

grain=grains(grains.findByLocation(g));
ori=grain.meanOrientation;


%transform line formed by the slip trace into crystal coordinate system
line=Miller(ori*line,grain.CS);

%Get image plane normal in crystal coordinates
implaneNormal=Miller(ori*zvector,grain.CS,'hkl');


%Define possible slip planes and slip directions for FCC structure
sPlane=symmetrise(Miller(1,1,1,grain.CS,'hkl'));

burgers=symmetrise(Miller(1,1,0,grain.CS,'uvw'),'antipodal');

%Calculate all possible slip trace vectors in crystal coordinates based on
%slip planes and SEM image plane normal
sTrace=cross(sPlane,implaneNormal);

%Get angle between the user-defined slip trace and the theoretically
%calculated slip trace directions
temp=angle(sTrace,line,'nosymmetry');

%Make lists of all possible slip planes and all trace-mismatch angles
allPlanes=[sPlane];
angles=[temp];


%Get the n slip planes that produce traces closest to the trace clicked on
%by the user.
[ASorted AIdx] = sort(angles);
if n<length(AIdx)
    top = AIdx(1:n);
else
    top=AIdx;
end


%Create list of the 3 most likely slip planes
activePlanes=allPlanes(top);

%Loop through list of 3 most likely slip planes.  For each one, calculate
%the possible slip directions.  Save the slip direction that has the
%largest component pointing out of the SEM image plane (since slip
%directions that lay in the image plane can't produce visible slip traces).
for i=1:length(top)
    [r,c] = find(isnull(dot_outer(vector3d(activePlanes(i)),vector3d(burgers))));
    activeBurgers=burgers(c);
    outOfPlaneAngle=abs(angle(activeBurgers,implaneNormal,'nosymmetry')-pi/2);
    activeBurgers=unique(activeBurgers(outOfPlaneAngle==min(outOfPlaneAngle)));

    plane(i)=activePlanes(i);
    direction(i)=activeBurgers;
    traceMismatch(i)=angles(top(i))/degree;
    slipDirAngle(i)=min(outOfPlaneAngle)/degree;
end


end

