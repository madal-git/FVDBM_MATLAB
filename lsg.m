function gradient = lsg(center_coor, center_value, sate_coor, sate_value)

%%%% gradient = lsg(center_coor, center_value, sate_coor, sate_value)
%%%% calculates the gradient at a given 2D location with a least-square
%%%% fashion.
%%%% center_coor is the 2D coordinates of the point where the gradient is
%%%% desired. A column vector
%%%% center_value is the value of the point where the gradient is desired.
%%%% sate_coor is the coordinate matrix of the known point. a 2 by n
%%%% matrix.
%%%% sate_value is the value vector or matrix of the known points. A m by n
%%%% matrix where m could be 1

[d1,l1]=size(sate_coor);
[d2,l2]=size(sate_value);

if d1~=2
    error('The dimension of the coordinates must be two!');
end

if l1~=l2
    error('The coordinates of satellite points must have the same length as their values!');
end

%% Calculate the coefficients
distance_vector_matrix=sate_coor-[ones(1,l1)*center_coor(1,1);ones(1,l1)*center_coor(2,1)];
distance_vector_matrix_sq=distance_vector_matrix.*distance_vector_matrix;
distance_sq=distance_vector_matrix_sq(1,:)+distance_vector_matrix_sq(2,:);
w=1./distance_sq;

%% Calculate the denomenator
 A=zeros(d1,d1);
 for j=1:l1
     dis_vec=sate_coor(:,j)-center_coor;
     A=A+w(j)*(dis_vec*dis_vec');
 end
%% Calculate the numerator
B=zeros(d1,d2);
for i=1:d2
    for j=1:l1
        dis_vec=sate_coor(:,j)-center_coor;
        B(:,i)=B(:,i)+w(j)*dis_vec*(sate_value(i,j)-center_value(i,1));
    end
end

 
 %% Calculate the gradient at the center point
 gradient=zeros(d1,d2);
 for i=1:d2
     gradient(:,i)=A\B(:,i);
 end