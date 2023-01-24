%% CPS assesment without edge-weight cut
clc
clearvars
close all

nlist_craft = readtable("./sim_data/craft_4players8_nlist.csv");
elist_craft = readtable("./sim_data/craft_4players8_elist.csv");

nlist_assembly = readtable("./sim_data/assembly_line_4players8_nlist.csv");
elist_assembly = readtable("./sim_data/assembly_line_4players8_elist.csv");

cps_craft = graph_sim_cps(nlist_craft,elist_craft);
cps_aseembly = graph_sim_cps(nlist_assembly,elist_assembly);

t1 = unique(elist_craft.t);
t2 = unique(elist_assembly.t);
ta = cps_aseembly;

figure('Color','w','Position',[50,50,1000,540])
tle = tiledlayout(2,3);
xlabel(tle,"Time (s)", 'FontWeight','bold')

nexttile
plot(cps_craft(:,1,2),'LineWidth',1.5,'Color','red')
xline(ta)
hold on
plot(cps_aseembly(:,1,2),'LineWidth',1.5,'Color','blue')
xlim('tight')
ylim([0,150])
ylabel("Distance (P1&P2)",'FontWeight','bold')
grid on
set(gca,'FontWeight','bold','LineWidth',1.5)
hold off

nexttile
plot(cps_craft(:,2,2),'LineWidth',1.5,'Color','red')
hold on
plot(cps_aseembly(:,2,2),'LineWidth',1.5,'Color','blue')
xlim('tight')
ylim([0,150])
ylabel("Distance (P1&P3)",'FontWeight','bold')
grid on
set(gca,'FontWeight','bold','LineWidth',1.5)
hold off

nexttile
plot(cps_craft(:,3,2),'LineWidth',1.5,'Color','red')
hold on
plot(cps_aseembly(:,3,2),'LineWidth',1.5,'Color','blue')
xlim('tight')
ylim([0,150])
ylabel("Distance (P1&P4)",'FontWeight','bold')
grid on
set(gca,'FontWeight','bold','LineWidth',1.5)
hold off

nexttile
plot(cps_craft(:,4,2),'LineWidth',1.5,'Color','red')
hold on
plot(cps_aseembly(:,4,2),'LineWidth',1.5,'Color','blue')
xlim('tight')
ylim([0,150])
ylabel("Distance (P2&P3)",'FontWeight','bold')
grid on
set(gca,'FontWeight','bold','LineWidth',1.5)
hold off

nexttile
plot(cps_craft(:,5,2),'LineWidth',1.5,'Color','red')
hold on
plot(cps_aseembly(:,5,2),'LineWidth',1.5,'Color','blue')
xlim('tight')
ylim([0,150])
ylabel("Distance (P2&P4)",'FontWeight','bold')
grid on
set(gca,'FontWeight','bold','LineWidth',1.5)
hold off

nexttile
plot(cps_craft(:,6,2),'LineWidth',1.5,'Color','red')
hold on
plot(cps_aseembly(:,6,2),'LineWidth',1.5,'Color','blue')
xlim('tight')
ylim([0,150])
ylabel("Distance (P3&P4)",'FontWeight','bold')
grid on
set(gca,'FontWeight','bold','LineWidth',1.5)
hold off

%% CPS assesment with edge-weight cut
clc
clearvars
close all

nlist = readtable("./sim_data/craft_4players8_nlist.csv");
elist = readtable("./sim_data/craft_4players8_elist.csv");

nlist.Properties.VariableNames = {'Name','Label'};
nlist.Name = string(nlist.Name);
num_nodes = size(nlist,1);
num_edges = size(elist,1);

edge_wcutoff = 16;
num_players = (size(elist,2)-1)/2;

% initialize graphs
G = cell(4,1);
A = cell(4,1);
for i=1:num_players
    G{i} = digraph(zeros(num_nodes,num_nodes),nlist);
    A{i} = adjacency(G{i});
end

out_mx = NaN(num_edges,nchoosek(num_players,2),2); % num_edges x 4C2 x 2 (adj_spectl,feat_sim)
for i=1:num_edges
    r = table2array(elist(i,:));
    gfeat = NaN(num_players,num_nodes,7);

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
        % adjacency matrices for spectral distance
        A{j} = adjacency(G{j});

        % graph features for feature-based distance        
        gfeat(j,:,1) = centrality(G{j},'outdegree','Importance',G{j}.Edges.Weight);
        gfeat(j,:,2) = centrality(G{j},'indegree','Importance',G{j}.Edges.Weight);
        gfeat(j,:,3) = centrality(G{j},'outcloseness');
        gfeat(j,:,4) = centrality(G{j},'incloseness');
        gfeat(j,:,5) = centrality(G{j},'betweenness');
        gfeat(j,:,6) = centrality(G{j},'pagerank','Importance',G{j}.Edges.Weight);
        gfeat(j,:,7) = centrality(G{j},'hubs','Importance',G{j}.Edges.Weight);
    end
    
    % graph feature & similarity
    l = 0;
    for j=1:(num_players-1)
        lambdaj = eigs(A{j});
        n = 1+l;
        for k=(j+1):(num_players)
            lambdak = eigs(A{k});
            out_mx(i,n,1) = norm(lambdaj-lambdak,2); % spectral
            out_mx(i,n,2) = norm(gfeat(j,:,:)-gfeat(k,:,:),'fro'); % feature based
            n = n+1;
        end
        l = l+(num_players-j);
    end
end

figure('Color','w','Position',[0,0,1000,540])
tle = tiledlayout(2,3);
xlabel(tle,"Time (s)", 'FontWeight','bold')

nexttile
plot(out_mx(:,1,2),'LineWidth',1.5,'Color','black')
ylim([0,60])
ylabel("Distance (P1&P2)",'FontWeight','bold')

nexttile
plot(out_mx(:,2,2),'LineWidth',1.5,'Color','black')
ylim([0,60])
ylabel("Distance (P1&P3)",'FontWeight','bold')

nexttile
plot(out_mx(:,3,2),'LineWidth',1.5,'Color','black')
ylim([0,60])
ylabel("Distance (P1&P4)",'FontWeight','bold')

nexttile
plot(out_mx(:,4,2),'LineWidth',1.5,'Color','black')
ylim([0,60])
ylabel("Distance (P2&P3)",'FontWeight','bold')

nexttile
plot(out_mx(:,5,2),'LineWidth',1.5,'Color','black')
ylim([0,60])
ylabel("Distance (P2&P4)",'FontWeight','bold')

nexttile
plot(out_mx(:,6,2),'LineWidth',1.5,'Color','black')
ylim([0,60])
ylabel("Distance (P3&P4)",'FontWeight','bold')

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

%% helper functions
function out_mx = graph_sim_cps(nlist,elist)
    tlist = unique(elist.t);
    nlist.Properties.VariableNames = {'Name','Label'};
    nlist.Name = string(nlist.Name);
    num_nodes = size(nlist,1);
    num_edges = length(tlist);
    num_players = (size(elist,2)-1)/2;
    
    % initialize graphs
    G = cell(4,1);
    A = cell(4,1);
    for i=1:num_players
        G{i} = digraph(zeros(num_nodes,num_nodes),nlist);
        A{i} = adjacency(G{i});
    end
    
    out_mx = NaN(num_edges,nchoosek(num_players,2),2); % num_edges x 4C2 x 2 (adj_spectl,feat_sim)
    for i=1:num_edges
        runningt = tlist(i);
        r = table2array(elist(elist.t==runningt,:));
        gfeat = NaN(num_players,num_nodes,7);
        
        % update graphs
        for j=1:size(r,1)
            for k=1:num_players
                s = r(j,2*k);
                t = r(j,2*k+1);
        
                if (s==0)&&(t==0)
                    % no update 
                elseif xor(s==0,t==0)
                    error("Invalid edge")
                else
                    [edge_exist,update_idx] = ismember([s,t],str2double(G{k}.Edges.EndNodes),'rows');
                    if (all(edge_exist))
                        w = G{k}.Edges.Weight;
                        G{k}.Edges.Weight(update_idx) = w(update_idx)+1;
                    else
                        G{k} = addedge(G{k},num2str(s),num2str(t),1);
                    end
                end
            end
        end
        
        % calculate graph features and CPS measure
        for j=1:num_players
            % adjacency matrices for spectral distance
            A{j} = adjacency(G{j});
    
            % graph features for feature-based distance        
            gfeat(j,:,1) = centrality(G{j},'outdegree','Importance',G{j}.Edges.Weight);
            gfeat(j,:,2) = centrality(G{j},'indegree','Importance',G{j}.Edges.Weight);
            gfeat(j,:,3) = centrality(G{j},'outcloseness');
            gfeat(j,:,4) = centrality(G{j},'incloseness');
            gfeat(j,:,5) = centrality(G{j},'betweenness');
            gfeat(j,:,6) = centrality(G{j},'pagerank','Importance',G{j}.Edges.Weight);
            gfeat(j,:,7) = centrality(G{j},'hubs','Importance',G{j}.Edges.Weight);
            
            % graph feature & similarity
            l = 0;
            for k=1:(num_players-1)
                lambdak = eigs(A{k});
                n = 1+l;
                for m=(k+1):(num_players)
                    lambdam = eigs(A{m});
                    out_mx(i,n,1) = norm(lambdak-lambdam,2); % spectral
                    out_mx(i,n,2) = norm(gfeat(k,:,:)-gfeat(m,:,:),'fro'); % feature based
                    n = n+1;
                end
                l = l+(num_players-k);
            end
        end
    end
end
