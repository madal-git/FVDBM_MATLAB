%% Change the cell info for the cell outside of the outer radius
for i=1:M
    CL=CELL{i};
    centriod_coord=CL{5};
    if dis(centriod_coord,[1;1])>(2/3)
        CL{10}=1;
        CL{11}=1;
        CL{12}=1;
        CL{19}=1;
        CL{20}=1;
        CL{21}=1;
        CL{37}=2;
    end
    CELL{i}=CL;
end

%% Change the node info for the cell outside of the outer radius
for i=1:N
    ND=NODE{i};
    cell_star=ND{5};
    for j=1:ND{4}
        CL=CELL{cell_star(j)};
        if CL{37}==2
            ND{2}=1;
            break;
        end
    end
    NODE{i}=ND;
end


%% Change the face info for the cell outside of the outer radius
for i=1:O
    FC=FACE{i};
    ND1=NODE{FC{8}};
    ND2=NODE{FC{9}};
    if ND1{2}==1 && ND2{2}==1
        FC{2}=1;
        FC{10}=1;
        FC{11}=1;
    end
    FACE{i}=FC;
end

%% Plot
figure
for r=1:M
    CL=CELL{r};
    centriod_coord=CL{5};
    plot(CL{22},CL{23},'black',CL{24},CL{25},'black',CL{26},CL{27},'black','linewidth',0.5);
    hold on
    if CL{37}==2
        plot(centriod_coord(1),centriod_coord(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'r','MarkerEdgeColor', 'r');
        hold on
    end
end

for r=1:N
    ND=NODE{r};
    nd_coord=ND{3};
    if ND{2}==1
        plot(nd_coord(1),nd_coord(2),'Marker', 'o','Markersize',2, 'MarkerFaceColor', 'g','MarkerEdgeColor', 'g');
        hold on
    end
end

for r=1:O
    FC=FACE{r};
    nd1_coord=FC{5};
    nd2_coord=FC{6};
    if FC{2}==1
        plot([nd1_coord(1),nd2_coord(1)],[nd1_coord(2),nd2_coord(2)],'blue','linewidth',1);
        hold on
    end
end

axis equal tight