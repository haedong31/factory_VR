%% Learning curves
clc
clearvars
close all

learning_curve = @(a,b,x) a.*power(x,-b);

r = 0.9:-0.1:0.5;
a = 1;
b = -(log(r)./log(2));
x = 1:100;

y = NaN(length(x),length(r));
for i=1:length(r)
    y(:,i) = learning_curve(a,b(i),x);
end

figure('Color','w')
newcolors = {'#F00','#F80','#FF0','#0B0','#00F'};
colororder(newcolors)
plot(x,y(:,1),'LineWidth',1.5)
hold on
for i=2:length(r)
    plot(x,y(:,i),'LineWidth',1.5)
end
hold off
xlabel("n^{th} unit to produce")
ylabel("Time")
legend(["90%","80%","70%","60%","50%"])
set(gca,'FontWeight','bold','LineWidth',1.5)

%% Snpashots of graphs
clc
clearvars
close all

exp_num = "exp1";
data_dir = fullfile("sim_data",exp_num);

nlist_craft = readtable(fullfile(data_dir,"craft8_nlist.csv"));
elist_craft = readtable(fullfile(data_dir,"craft8_elist.csv"));
nlist_mass = readtable(fullfile(data_dir,"mass8_nlist.csv"));
elist_mass = readtable(fullfile(data_dir,"mass8_elist.csv"));

nlist_craft.Properties.VariableNames = {'Name','Label'};
nlist_craft.Name = string(nlist_craft.Name);
nlist_mass.Properties.VariableNames = {'Name','Label'};
nlist_mass.Name = string(nlist_mass.Name);

% snap_pts = [floor(quantile(1:num_edges,0.25)),floor(quantile(1:num_edges,0.5)),...
%     floor(quantile(1:num_edges,0.75)),floor(quantile(1:num_edges,1))];

snap_pts = NaN(4,3);
snap_pts([1,3],:) = repmat([10,40,80],2,1);
snap_pts([2,4],:) = repmat([10,45,90],2,1);

edge_cutoff = 16;
Gsnaps_craft = take_snaps(nlist_craft,elist_craft,snap_pts,edge_cutoff);
Gsnaps_mass = take_snaps(nlist_mass,elist_mass,snap_pts,edge_cutoff);

f = figure('Color','w','Position',[50,50,1000,540]);
orient(f,'landscape')
tle = tiledlayout(2,3);
title(tle,"Snapshots P1&P3",'FontWeight','bold')
for i=1:3
    subplot(2,3,i)
    plot(Gsnaps_mass{1,i},'NodeLabel',nlist_mass.Label,'NodeColor','b','MarkerSize',10,...
        'NodeFontSize',11,'NodeFontWeight','bold',...
        'EdgeColor','b','LineWidth',Gsnaps_mass{1,i}.Edges.Weight,'Layout','auto')
    set(gca,'Visible','off')
    
    subplot(2,3,i+3)
    plot(Gsnaps_mass{3,i},'NodeLabel',nlist_mass.Label,'NodeColor','r','MarkerSize',10,...
        'NodeFontSize',11,'NodeFontWeight','bold',...
        'EdgeColor','r','LineWidth',Gsnaps_mass{3,i}.Edges.Weight,'Layout','auto')
    set(gca,'Visible','off')
end


%% CPS assesment
clc
clearvars
close all

exp_num = "exp1";
data_dir = fullfile("sim_data",exp_num);

minfo = readtable(fullfile(data_dir,"mass8_minfo.csv"));
nlist_craft = readtable(fullfile(data_dir,"craft8_nlist.csv"));
elist_craft = readtable(fullfile(data_dir,"craft8_elist.csv"));
nlist_assembly = readtable(fullfile(data_dir,"mass8_nlist.csv"));
elist_assembly = readtable(fullfile(data_dir,"mass8_elist.csv"));

cps_craft = graph_sim_cps(nlist_craft,elist_craft);
cps_aseembly = graph_sim_cps(nlist_assembly,elist_assembly);

t1 = unique(elist_craft.t);
t2 = unique(elist_assembly.t);
tearly = min(t1(end),t2(end));
tlate = max(t1(end),t2(end));

c = parula(4);

f = figure('Color','w','Position',[50,50,800,430]);
orient(f,'landscape')
tle = tiledlayout(2,3);
% title(tle,"Same Learning Curves (70% & 70%)",'FontWeight','bold')
xlabel(tle,"Time (m)", 'FontWeight','bold')

nexttile
sub_minfo = minfo(string(minfo.player)=="player1"&string(minfo.item)=="finish1",:);
specialt1 = sub_minfo.time_stamp;
specialt_idx1 = find_specialt_idx(specialt1,t2);

sub_minfo = minfo(string(minfo.player)=="player2"&string(minfo.item)=="base",:);
specialt2 = sub_minfo.time_stamp;
specialt_idx2 = find_specialt_idx(specialt2,t2);

sub_minfo = minfo(string(minfo.player)=="player2"&string(minfo.item)=="finish2",:);
specialt3 = sub_minfo.time_stamp;
specialt_idx3 = find_specialt_idx(specialt3,t2);

plot(t1,cps_craft(:,1,2),'--','LineWidth',1.5,'Color','red')
hold on
plot(t2,cps_aseembly(:,1,2),'LineWidth',1.5,'Color','blue')
scatter(t2(specialt_idx1),cps_aseembly(specialt_idx1,1,2),30,'filled', ...
    'MarkerFaceColor',c(1,:))
scatter(t2(specialt_idx2),cps_aseembly(specialt_idx2,1,2),30,'filled', ...
    'MarkerFaceColor',c(2,:))
scatter(t2(specialt_idx3),cps_aseembly(specialt_idx3,1,2),30,'filled', ...
    'MarkerFaceColor',c(3,:))
hold off
xline(tearly,'LineWidth',1.5,'Color','blue')
xline(tlate,'LineWidth',1.5,'Color','red')
xlim('tight')
ylim([0,150])
ylabel("CPS (P1&P2)",'FontWeight','bold')
grid on
legend(["Craft","Mass","P1 Module","P2 Submodule","P2 Module"],'Location','northwest')
set(gca,'FontWeight','bold','LineWidth',1.5)

% plot 2
nexttile
sub_minfo = minfo(string(minfo.player)=="player1"&string(minfo.item)=="finish1",:);
specialt1 = sub_minfo.time_stamp;
specialt_idx1 = find_specialt_idx(specialt1,t2);

sub_minfo = minfo(string(minfo.player)=="player3"&string(minfo.item)=="steer",:);
specialt2 = sub_minfo.time_stamp;
specialt_idx2 = find_specialt_idx(specialt2,t2);

sub_minfo = minfo(string(minfo.player)=="player3"&string(minfo.item)=="finish3",:);
specialt3 = sub_minfo.time_stamp;
specialt_idx3 = find_specialt_idx(specialt3,t2);

plot(t1,cps_craft(:,2,2),'--','LineWidth',1.5,'Color','red')
hold on
plot(t2,cps_aseembly(:,2,2),'LineWidth',1.5,'Color','blue')
scatter(t2(specialt_idx1),cps_aseembly(specialt_idx1,2,2),30,'filled', ...
    'MarkerFaceColor',c(1,:))
scatter(t2(specialt_idx2),cps_aseembly(specialt_idx2,2,2),30,'filled', ...
    'MarkerFaceColor',c(2,:))
scatter(t2(specialt_idx3),cps_aseembly(specialt_idx3,2,2),30,'filled', ...
    'MarkerFaceColor',c(3,:))
hold off
xlim('tight')
ylim([0,150])
ylabel("CPS (P1&P3)",'FontWeight','bold')
grid on
legend(["","","P1 Module","P3 Submodule","P3 Module"],'Location','northwest')
set(gca,'FontWeight','bold','LineWidth',1.5)

% plot 3
nexttile
sub_minfo = minfo(string(minfo.player)=="player1"&string(minfo.item)=="finish1",:);
specialt1 = sub_minfo.time_stamp;
specialt_idx1 = find_specialt_idx(specialt1,t2);

sub_minfo = minfo(string(minfo.player)=="player4"&string(minfo.item)=="body",:);
specialt2 = sub_minfo.time_stamp;
specialt_idx2 = find_specialt_idx(specialt2,t2);

sub_minfo = minfo(string(minfo.player)=="player4"&string(minfo.item)=="finish4",:);
specialt3 = sub_minfo.time_stamp;
specialt_idx3 = find_specialt_idx(specialt3,t2);

plot(t1,cps_craft(:,3,2),'--','LineWidth',1.5,'Color','red')
hold on
plot(t2,cps_aseembly(:,3,2),'LineWidth',1.5,'Color','blue')
scatter(t2(specialt_idx1),cps_aseembly(specialt_idx1,3,2),30,'filled', ...
    'MarkerFaceColor',c(1,:))
scatter(t2(specialt_idx2),cps_aseembly(specialt_idx2,3,2),30,'filled', ...
    'MarkerFaceColor',c(2,:))
scatter(t2(specialt_idx3),cps_aseembly(specialt_idx3,3,2),30,'filled', ...
    'MarkerFaceColor',c(3,:))
hold off
xlim('tight')
ylim([0,150])
ylabel("CPS (P1&P4)",'FontWeight','bold')
grid on
legend(["","","P1 Module","P4 Submodule","P4 Module"],'Location','northwest')
set(gca,'FontWeight','bold','LineWidth',1.5)

% plot 4
nexttile
sub_minfo = minfo(string(minfo.player)=="player2"&string(minfo.item)=="base",:);
specialt1 = sub_minfo.time_stamp;
specialt_idx1 = find_specialt_idx(specialt1,t2);

sub_minfo = minfo(string(minfo.player)=="player2"&string(minfo.item)=="finish2",:);
specialt2 = sub_minfo.time_stamp;
specialt_idx2 = find_specialt_idx(specialt2,t2);

sub_minfo = minfo(string(minfo.player)=="player3"&string(minfo.item)=="steer",:);
specialt3 = sub_minfo.time_stamp;
specialt_idx3 = find_specialt_idx(specialt3,t2);

sub_minfo = minfo(string(minfo.player)=="player3"&string(minfo.item)=="finish3",:);
specialt4 = sub_minfo.time_stamp;
specialt_idx4 = find_specialt_idx(specialt4,t2);

plot(t1,cps_craft(:,4,2),'--','LineWidth',1.5,'Color','red')
hold on
plot(t2,cps_aseembly(:,4,2),'LineWidth',1.5,'Color','blue')
scatter(t2(specialt_idx1),cps_aseembly(specialt_idx1,4,2),30,'filled', ...
    'MarkerFaceColor',c(1,:))
scatter(t2(specialt_idx2),cps_aseembly(specialt_idx2,4,2),30,'filled', ...
    'MarkerFaceColor',c(2,:))
scatter(t2(specialt_idx3),cps_aseembly(specialt_idx3,4,2),30,'filled', ...
    'MarkerFaceColor',c(3,:))
scatter(t2(specialt_idx4),cps_aseembly(specialt_idx4,4,2),30,'filled', ...
    'MarkerFaceColor',c(4,:))
hold off
xlim('tight')
ylim([0,150])
ylabel("CPS (P2&P3)",'FontWeight','bold')
grid on
legend(["","","P2 Submodule","P2 Module","P3 Submodule","P3 Module"],'Location','northwest')
set(gca,'FontWeight','bold','LineWidth',1.5)

% plot 5
nexttile
sub_minfo = minfo(string(minfo.player)=="player2"&string(minfo.item)=="base",:);
specialt1 = sub_minfo.time_stamp;
specialt_idx1 = find_specialt_idx(specialt1,t2);

sub_minfo = minfo(string(minfo.player)=="player2"&string(minfo.item)=="finish2",:);
specialt2 = sub_minfo.time_stamp;
specialt_idx2 = find_specialt_idx(specialt2,t2);

sub_minfo = minfo(string(minfo.player)=="player4"&string(minfo.item)=="body",:);
specialt3 = sub_minfo.time_stamp;
specialt_idx3 = find_specialt_idx(specialt3,t2);

sub_minfo = minfo(string(minfo.player)=="player4"&string(minfo.item)=="finish4",:);
specialt4 = sub_minfo.time_stamp;
specialt_idx4 = find_specialt_idx(specialt4,t2);

plot(t1,cps_craft(:,5,2),'--','LineWidth',1.5,'Color','red')
hold on
plot(t2,cps_aseembly(:,5,2),'LineWidth',1.5,'Color','blue')
scatter(t2(specialt_idx1),cps_aseembly(specialt_idx1,5,2),30,'filled', ...
    'MarkerFaceColor',c(1,:))
scatter(t2(specialt_idx2),cps_aseembly(specialt_idx2,5,2),30,'filled', ...
    'MarkerFaceColor',c(2,:))
scatter(t2(specialt_idx3),cps_aseembly(specialt_idx3,5,2),30,'filled', ...
    'MarkerFaceColor',c(3,:))
scatter(t2(specialt_idx4),cps_aseembly(specialt_idx4,5,2),30,'filled', ...
    'MarkerFaceColor',c(4,:))
xlim('tight')
ylim([0,150])
ylabel("CPS (P2&P4)",'FontWeight','bold')
grid on
legend(["","","P2 Submodule","P2 Module","P4 Submodule","P4 Module"],'Location','northwest')
set(gca,'FontWeight','bold','LineWidth',1.5)
hold off

% plot 6
nexttile
sub_minfo = minfo(string(minfo.player)=="player3"&string(minfo.item)=="steer",:);
specialt1 = sub_minfo.time_stamp;
specialt_idx1 = find_specialt_idx(specialt1,t2);

sub_minfo = minfo(string(minfo.player)=="player3"&string(minfo.item)=="finish3",:);
specialt2 = sub_minfo.time_stamp;
specialt_idx2 = find_specialt_idx(specialt2,t2);

sub_minfo = minfo(string(minfo.player)=="player4"&string(minfo.item)=="body",:);
specialt3 = sub_minfo.time_stamp;
specialt_idx3 = find_specialt_idx(specialt3,t2);

sub_minfo = minfo(string(minfo.player)=="player4"&string(minfo.item)=="finish4",:);
specialt4 = sub_minfo.time_stamp;
specialt_idx4 = find_specialt_idx(specialt4,t2);

plot(t1,cps_craft(:,6,2),'--','LineWidth',1.5,'Color','red')
hold on
plot(t2,cps_aseembly(:,6,2),'LineWidth',1.5,'Color','blue')
scatter(t2(specialt_idx1),cps_aseembly(specialt_idx1,6,2),30,'filled', ...
    'MarkerFaceColor',c(1,:))
scatter(t2(specialt_idx2),cps_aseembly(specialt_idx2,6,2),30,'filled', ...
    'MarkerFaceColor',c(2,:))
scatter(t2(specialt_idx3),cps_aseembly(specialt_idx3,6,2),30,'filled', ...
    'MarkerFaceColor',c(3,:))
scatter(t2(specialt_idx4),cps_aseembly(specialt_idx4,6,2),30,'filled', ...
    'MarkerFaceColor',c(4,:))
xlim('tight')
ylim([0,150])
ylabel("CPS (P3&P4)",'FontWeight','bold')
grid on
legend(["","","P3 Submodule","P3 Module","P4 Submodule","P4 Module"],'Location','northwest')
set(gca,'FontWeight','bold','LineWidth',1.5)
hold off

%% Low-dimensional embedding for user-action space
clc
close all
clearvars

nlist_craft = readtable("./sim_data/craft_4players8_nlist.csv");
elist_craft = readtable("./sim_data/craft_4players8_elist.csv");
nlist_assembly = readtable("./sim_data/assembly_line_4players8_nlist.csv");
elist_assembly = readtable("./sim_data/assembly_line_4players8_elist.csv");

gfeat_craft = graph_feat(nlist_craft,elist_craft);
gfeat_assembly = graph_feat(nlist_assembly,elist_assembly);
acts = ["Spawn","Move","Assemble","Pass"];

%----- craft production -----%
dims = size(gfeat_craft);
gfeat_craft2 = reshape(gfeat_craft,[dims(1)*dims(2)*dims(3),dims(4)]);

% % label 1: category-level 1
% l1 = repmat(unique(elist_craft.t),[dims(2)*dims(3),1]);
% 
% % label 2: category-level 2
% l2_unit = strings(dims(1)*dims(2),1);
% for i=1:dims(2)
%     l2_unit(1+(i-1)*dims(1):i*dims(1)) = strcat("Player",num2str(i));
% end
% l2 = strings(dims(1)*dims(2)*dims(3),1);
% for i=1:dims(3)
%     l2(1+(i-1)*dims(1)*dims(2):i*dims(1)*dims(2)) = l2_unit;
% end

% label 3: category-level 3
l3 = strings(dims(1)*dims(2)*dims(3),1);
for i=1:dims(3)
    l3(1+(i-1)*dims(1)*dims(2):i*dims(1)*dims(2)) = acts(i);
end

c = NaN(length(l3),3);
idx1 = (l3=="Spawn");
idx2 = (l3=="Move");
idx3 = (l3=="Assemble");
idx4 = (l3=="Pass");

cmap_end = parula(4);
cmap_begin = [0.7, 0.65, 1;
    0.7, 0.87, 1;
    0.7, 0.8, 0.6;
    0.98, 0.99, 0.7];
cmap = [linspace(cmap_begin(1,1),cmap_end(1,1),dims(1))', ...
    linspace(cmap_begin(1,2),cmap_end(1,2),dims(1))', ...
    linspace(cmap_begin(1,3),cmap_end(1,3),dims(1))'];
c(idx1,:) = repmat(cmap,dims(2),1);
cmap = [linspace(cmap_begin(2,1),cmap_end(2,1),dims(1))', ...
    linspace(cmap_begin(2,2),cmap_end(2,2),dims(1))', ...
    linspace(cmap_begin(2,3),cmap_end(2,3),dims(1))'];
c(idx2,:) = repmat(cmap,dims(2),1);
cmap = [linspace(cmap_begin(3,1),cmap_end(3,1),dims(1))', ...
    linspace(cmap_begin(3,2),cmap_end(3,2),dims(1))', ...
    linspace(cmap_begin(3,3),cmap_end(3,3),dims(1))'];
c(idx3,:) = repmat(cmap,dims(2),1);
cmap = [linspace(cmap_begin(4,1),cmap_end(4,1),dims(1))', ...
    linspace(cmap_begin(4,2),cmap_end(4,2),dims(1))', ...
    linspace(cmap_begin(4,3),cmap_end(4,3),dims(1))'];
c(idx4,:) = repmat(cmap,dims(2),1);

rng(231)
ops = statset('MaxIter',2000);
f = figure('Color','w','Position',[100,100,920,320]);
orient(f,'landscape')

subplot(1,2,1)
embd = tsne(gfeat_craft2,'Distance','jaccard','NumDimensions',2, ...
    'Standardize',true,'Options',ops);
scatter(embd(idx1,1),embd(idx1,2),'o','filled','CData',c(idx1,:))
hold on
scatter(embd(idx2,1),embd(idx2,2),'square','filled','CData',c(idx2,:))
scatter(embd(idx3,1),embd(idx3,2),'diamond','filled','CData',c(idx3,:))
scatter(embd(idx4,1),embd(idx4,2),'^','filled','CData',c(idx4,:))
hold off
xlabel("Embedding Dim1")
ylabel("Embedding Dim2")
axis tight
legend(acts)
set(gca,'FontWeight','bold','LineWidth',1.5)

%----- mass production -----%
dims = size(gfeat_assembly);
gfeat_assembly2 = reshape(gfeat_assembly,[dims(1)*dims(2)*dims(3),dims(4)]);

l3 = strings(dims(1)*dims(2)*dims(3),1);
for i=1:dims(3)
    l3(1+(i-1)*dims(1)*dims(2):i*dims(1)*dims(2)) = acts(i);
end

c = NaN(length(l3),3);
idx1 = (l3=="Spawn");
idx2 = (l3=="Move");
idx3 = (l3=="Assemble");
idx4 = (l3=="Pass");

cmap = [linspace(cmap_begin(1,1),cmap_end(1,1),dims(1))', ...
    linspace(cmap_begin(1,2),cmap_end(1,2),dims(1))', ...
    linspace(cmap_begin(1,3),cmap_end(1,3),dims(1))'];
c(idx1,:) = repmat(cmap,dims(2),1);
cmap = [linspace(cmap_begin(2,1),cmap_end(2,1),dims(1))', ...
    linspace(cmap_begin(2,2),cmap_end(2,2),dims(1))', ...
    linspace(cmap_begin(2,3),cmap_end(2,3),dims(1))'];
c(idx2,:) = repmat(cmap,dims(2),1);
cmap = [linspace(cmap_begin(3,1),cmap_end(3,1),dims(1))', ...
    linspace(cmap_begin(3,2),cmap_end(3,2),dims(1))', ...
    linspace(cmap_begin(3,3),cmap_end(3,3),dims(1))'];
c(idx3,:) = repmat(cmap,dims(2),1);
cmap = [linspace(cmap_begin(4,1),cmap_end(4,1),dims(1))', ...
    linspace(cmap_begin(4,2),cmap_end(4,2),dims(1))', ...
    linspace(cmap_begin(4,3),cmap_end(4,3),dims(1))'];
c(idx4,:) = repmat(cmap,dims(2),1);

subplot(1,2,2)
embd = tsne(gfeat_assembly2,'Distance','jaccard','NumDimensions',2, ...
    'Standardize',true,'Options',ops);
scatter(embd(idx1,1),embd(idx1,2),'o','filled','CData',c(idx1,:))
hold on
scatter(embd(idx2,1),embd(idx2,2),'square','filled','CData',c(idx2,:))
scatter(embd(idx3,1),embd(idx3,2),'diamond','filled','CData',c(idx3,:))
scatter(embd(idx4,1),embd(idx4,2),'^','filled','CData',c(idx4,:))
hold off
xlabel("Embedding Dim1")
ylabel("Embedding Dim2")
axis tight
set(gca,'FontWeight','bold','LineWidth',1.5)

%% helper functions
function out_mx = graph_sim_cps(nlist,elist)
    tlist = unique(elist.t);
    nlist.Properties.VariableNames = {'Name','Label'};
    nlist.Name = string(nlist.Name);
    num_nodes = size(nlist,1);
    num_tframes = length(tlist);
    num_players = (size(elist,2)-1)/2;
    
    % initialize graphs
    G = cell(4,1);
    A = cell(4,1);
    for i=1:num_players
        G{i} = digraph(zeros(num_nodes,num_nodes),nlist);
        A{i} = adjacency(G{i});
    end
    
    out_mx = NaN(num_tframes,nchoosek(num_players,2),2); % num_tframes x 4C2 x 2 (adj_spectl,feat_sim)
    for i=1:num_tframes
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
            gfeat(j,:,6) = centrality(G{j},'pagerank','Importance',G{j}.Edges.Weight,'MaxIterations',1000);
            gfeat(j,:,7) = centrality(G{j},'hubs','Importance',G{j}.Edges.Weight,'MaxIterations',1000);
            
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

function gfeat = graph_feat(nlist,elist)
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
    
    gfeat = NaN(num_edges,num_players,num_nodes,7); % t x p x a x f
    for i=1:num_edges
        runningt = tlist(i);
        r = table2array(elist(elist.t==runningt,:));
        
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
            gfeat(i,j,:,1) = centrality(G{j},'outdegree','Importance',G{j}.Edges.Weight);
            gfeat(i,j,:,2) = centrality(G{j},'indegree','Importance',G{j}.Edges.Weight);
            gfeat(i,j,:,3) = centrality(G{j},'outcloseness');
            gfeat(i,j,:,4) = centrality(G{j},'incloseness');
            gfeat(i,j,:,5) = centrality(G{j},'betweenness');
            gfeat(i,j,:,6) = centrality(G{j},'pagerank','Importance',G{j}.Edges.Weight);
            gfeat(i,j,:,7) = centrality(G{j},'hubs','Importance',G{j}.Edges.Weight);  
        end
    end
end

function specialt_idx = find_specialt_idx(t,tt)
    specialt_idx = NaN(length(t),1);
    for i=1:length(t)
        [~,specialt_idx(i)] = min(abs(tt-t(i)));
    end
end

function Gsnaps = take_snaps(nlist,elist,snap_pts,edge_wcutoff)
    tlist = unique(elist.t);
    nlist.Properties.VariableNames = {'Name','Label'};
    nlist.Name = string(nlist.Name);    

    num_nodes = size(nlist,1);
    num_tframes = length(tlist);
    [num_players,num_snaps] = size(snap_pts);
    
    G = cell(num_players,1);
    Gsnaps = cell(num_players,num_snaps);
    for i=1:num_players
        G{i} = digraph(zeros(num_nodes,num_nodes),nlist);
    end
    
    % have graphs evolve
    tpointer = ones(num_players,1);
    for i=1:num_tframes
        runningt = tlist(i);
        r = table2array(elist(elist.t==runningt,:));

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
                        if (w(update_idx) >= edge_wcutoff)
                            G{k}.Edges.Weight(update_idx) = w(update_idx)+(0.5)^(w(update_idx)-edge_wcutoff+1);
                        else
                            G{k}.Edges.Weight(update_idx) = w(update_idx)+1;
                        end
                    else
                        G{k} = addedge(G{k},num2str(s),num2str(t),1);
                    end
                end
            
                if tpointer(k) <= num_snaps
                    if i==snap_pts(k,tpointer(k))
                        Gsnaps{k,tpointer(k)} = G{k};
                        tpointer(k) = tpointer(k)+1;
                    end
                end
            end
        end
    end
end