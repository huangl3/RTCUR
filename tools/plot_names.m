function g = plot_names(tmp,author,cidx)
g1=find(cidx==3);
g2=find(cidx==2);
g3=find(cidx==1);

key_author=[%"David Dunson","Lawrence Carin","Yongtao Guan",...
    %"Brian J Reich","Song Xi Chen","Fang Yao","Ying Wei","Heping Zhang","Harrison H Zhou",...
%     "Debashis Paul","Yi Lin","Michael J Todd","Kani Chen","Qiwei Yao","Qi-Man Shao","Ming Yuan",...
%     "Richard Samworth","J S Marron","Ming-Hui Chen","Yanyuan Ma","Cheng Yong Tang","Bing-Yi Jing",...
%     "Michael Sherman","Rasmus Waagepetersen","Cun-Hui Zhang","Hongzhe Li","Jian Huang",...
%     "David L Donoho","Grace Wahba","Annie Qu","Lawrence D Brown","Rui Song","Xiaotong Shen",...
   "Yichao Wu","Guang Cheng","Trevor J Hastie","Zhezhen Jin","Wang Zhou",...
   "Kung-Sik Chan","Lan Zhang","Jacob Bien","Jiancheng Jiang","Daniela M Witten",...
    "Hansheng Wang","Yacine Ait-Sahalia","Noelle I Samia","Dan Yu Lin","Hsin-Cheng Huang",...
    "Hua Liang","Ming-Yen Cheng","Xiang Liu","Clifford Lam","Qiwei Yao","Lan Zhou","Richard Samworth"...
    "T Tony Cai", "Jiashun Jin", "Amnon Neeman","Haibo Zhou","Jiancheng Jiang","Yingying Fan","Tiejun Tong",...
    "Yi-Hau Chen","Philip J Brown","Jane-Ling Wang","Hans-Georg Muller","Yanyuan Ma","Yuedong Wang"...
    "Hui Zou","Chenlei Leng","Yingcun Xia","Qingxia Chen","Byeong U Park","Bo Li","Stefan Sperlich"...
    "Wei Pan","Joseph G Ibrahim","Heng Peng","Runze Li","Philip J Brown","Yi-Hau Chen","Nilanjan Chatterjee",...
    "Jianqing Fan","Raymond J Carroll","Peter Hall","Qingxia Chen","Joseph G Ibrahim","Ming-Hui Chen"% "Hua Liang","Tony Cai","Hans-Georg Muller","Jing Qin"...
    ];
toshow_id=find(contains(author,key_author));
pU=cmdscale(squareform(pdist(tmp)),2);
extra_id1 = find(pU(:,1)>10 & pU(:,1)<22);
extra_id2 = find(pU(:,2)>10);
extra_id3 = find(pU(:,1)<-12);
toshow_id = [toshow_id;extra_id1;extra_id2;extra_id3];
toshowg1=intersect(g1,toshow_id);
toshowg2=intersect(g2,toshow_id);
toshowg3=intersect(g3,toshow_id);

pU(g1,:) = pU(g1,:)+0.5*randn(size(pU(g1,:)));
pU(g2,:) = pU(g2,:)+0.5*randn(size(pU(g2,:)));
pU(g3,:) = pU(g3,:)+0.5*randn(size(pU(g3,:)));
figure;
g = plot(pU(g1,1),pU(g1,2),'.','color',[0, 0.4470, 0.7410],'MarkerSize',14);
hold on;
plot(pU(g2,1),pU(g2,2),'.','color',	[0.8500, 0.3250, 0.0980],'MarkerSize',14)
plot(pU(g3,1),pU(g3,2),'.','color',[0.4660, 0.6740, 0.1880],'MarkerSize',14)
para1.color = [0, 0.4470, 0.7410]; para1.fontsize = 14;
para2.color = [0.8500, 0.3250, 0.0980]; para2.fontsize = 14;
para3.color = [0.4660, 0.6740, 0.1880]; para3.fontsize = 14;
textfit(pU(toshowg1,1),pU(toshowg1,2),cellstr(author(toshowg1)),para1)
textfit(pU(toshowg2,1),pU(toshowg2,2),cellstr(author(toshowg2)),para2)
textfit(pU(toshowg3,1),pU(toshowg3,2),cellstr(author(toshowg3)),para3)
% set(gcf,'Position',[1000 593 900 750])
xlabel("First Principle Component of SCORE result",'fontsize',14);
ylabel("Second Principle Component of SCORE result",'fontsize',14);
set(gcf,'renderer','Painters')
end
