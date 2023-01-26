%% Dynamic graphs snapshots (4 times)
clc
clearvars
close all

nlist = readtable("./sim_data/assembly_line_4players8_nlist.csv");
elist = readtable("./sim_data/assembly_line_4players8_elist.csv");

nlist.Properties.VariableNames = {'Name','Label'};
nlist.Name = string(nlist.Name);
num_nodes = size(nlist,1);
num_edges = size(elist,1);

edge_wcutoff = 16;
num_players = (size(elist,2)-1)/2;

% initialize graphs; snap graphs at 4 time points
G = cell(4,1);
Gsnaps = cell(4,4);
for i=1:num_players
    G{i} = digraph(zeros(num_nodes,num_nodes),nlist);
end
snap_pts = [floor(quantile(1:num_edges,0.25)),floor(quantile(1:num_edges,0.5)),...
    floor(quantile(1:num_edges,0.75)),floor(quantile(1:num_edges,1))];

% have graphs evolve
tpointer = 1;
for i=1:num_edges
    r = table2array(elist(i,:));
    for j=1:num_players
        s = r(2*j);
        t = r(2*j+1);

        if (s==0)&&(t==0)
            % no update 
        elseif xor(s==0,t==0)
            error("Invalid edge")
        else
            [edge_exist,update_idx] = ismember([s,t],str2double(G{j}.Edges.EndNodes),'rows');
            if (all(edge_exist))
                w = G{j}.Edges.Weight;
                if (w(update_idx) >= edge_wcutoff)
                    G{j}.Edges.Weight(update_idx) = w(update_idx)+(0.5)^(w(update_idx)-edge_wcutoff+1);
                else
                    G{j}.Edges.Weight(update_idx) = w(update_idx)+1;
                end
            else
                G{j} = addedge(G{j},num2str(s),num2str(t),1);
            end
        end
    end

    if i==snap_pts(tpointer)
        for j=1:num_players
            Gsnaps{j,tpointer} = G{j};
        end
        tpointer = tpointer+1;
    end
end

ccodes = parula(num_players);
figure('Color','w','Position',[10,10,1200,1300])
for i=1:num_players
    for j=1:length(snap_pts)
        subplot(4,4,num_players*(i-1)+j)
        plot(Gsnaps{i,j},'NodeLabel',nlist.Label,'NodeColor',ccodes(i,:),'MarkerSize',10,...
            'NodeFontSize',11,'NodeFontWeight','bold',...
            'EdgeColor',ccodes(i,:),'LineWidth',Gsnaps{i,j}.Edges.Weight,'Layout','auto')
        set(gca,'Visible','off')
    end
end