all_pvals = rsid_pval
expected = np.random.uniform(0,1,all_pvals.shape[0])

bordercolor = '#333333'
borderwidth = 5
axis_font_size = 30
label_font_size = 25
legend_font_size = 25

fig = plt.figure(figsize=(20,20))
ax = fig.add_subplot(111)

for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(20)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(20)
    tick.label1.set_fontweight('bold')

#plt.xticks(axis_midpoint, axis_label, rotation='vertical')
ax.tick_params(axis='both', which = 'major',
                   length = 10, width = borderwidth, pad=10,
                   labelsize = label_font_size+10,
                   color = bordercolor,
                   labelcolor = bordercolor,
                   bottom = True, top = False, left = True, right = False)

               
xval = -np.log10(expected[np.argsort(expected)])
yval = -np.log10(all_pvals[np.argsort(all_pvals)])

plt.scatter(xval, yval, color="#C10020")
ax.set_title("P-P plot of the pvalues from expected and \n observed pvalues for association test", fontsize=30,fontweight="bold")
ax.plot([0,6],[0,6],'r--',color="black",lw=2)
ax.set_xlim(0,6)
ax.set_ylim(0,6)
for side, border in ax.spines.items():
    border.set_linewidth(borderwidth)
    border.set_color(bordercolor)
              
ax.set_xlabel(r"expected $-\log_{10}(Pval)$",{'size': axis_font_size+10, 'color': bordercolor},labelpad = 25, fontweight="bold")
ax.set_ylabel(r"observed $-\log_{10}(Pval)$",{'size': axis_font_size+10, 'color': bordercolor},labelpad = 25, fontweight="bold")
plt.savefig("P-P plot of the pvalues from expected and observed after modifying exponential")
plt.show()

#############################################  Manhattan PLot                        ########################################
col = ["wheat","silver"]
#pcol = ["#C10020","#803E75","#C10020","#803E75","#C10020","#803E75","#C10020","#803E75","#C10020","#803E75","#C10020","#803E75","#C10020","#803E75","#C10020","#803E75","#C10020","#803E75","#C10020","#803E75","#C10020","#803E75"]
pcol=["#C10020","blue","#C10020","blue","#C10020","blue","#C10020","blue","#C10020","blue","#C10020","blue","#C10020","blue","#C10020","blue","#C10020","blue","#C10020","blue","#C10020","blue","#C10020","blue","#C10020","blue","#C10020","blue","#C10020","blue","#C10020","blue","#C10020","blue","#C10020","blue","#C10020","blue","#C10020","blue","#C10020","blue","#C10020","blue"]
colored_indices = np.where(all_pvals>4.3)[0]
bp_c = all_basepairs[colored_indices]
rsid_c = all_pvals[colored_indices]

rang = list()
for i  in range(len(rsid_info)-1):
    rang.append(range(rsid_info[i],rsid_info[i+1]))
#print(rang)



bordercolor = '#333333'
borderwidth = 5
axis_font_size = 30
label_font_size = 25
legend_font_size = 25

fig = plt.figure(figsize=(40,20))
ax = fig.add_subplot(111)
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(20)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(20)
    tick.label1.set_fontweight('bold')
axis_midpoint = [(i/2)+np.sum(last_axis[0:j]) if j!=0 else i/2 for j,i in enumerate(last_axis)]
axis_label = ["chr{0}".format(i) for i in range(1,23)]
plt.xticks(axis_midpoint, axis_label, rotation='vertical')
ax.tick_params(axis='both', which = 'major',
                   length = 10, width = borderwidth, pad=10,
                   labelsize = label_font_size+10,
                   color = bordercolor,
                   labelcolor = bordercolor,
                   bottom = True, top = False, left = True, right = False
                  )
for i in range(len(rsid_info)):
    if i < 22:
        if i%2 == 0:
            colors=col[0]
            ax.scatter(all_basepairs[rsid_info[i]:rsid_info[i+1]], all_pvals[rsid_info[i]:rsid_info[i+1]], color=colors, s = 2, alpha=0.04, linewidths=4)
        else:
            colors=col[1]
            ax.scatter(all_basepairs[rsid_info[i]:rsid_info[i+1]], all_pvals[rsid_info[i]:rsid_info[i+1]], color=colors, s = 2, alpha=0.04, linewidths=4)
    else:
        colors=col[0]
        ax.scatter(all_basepairs[rsid_info[i]:], all_pvals[rsid_info[i]:], color=colors, s = 2, alpha=0.04, linewidths=4)

for i,j in enumerate(colored_indices):
    s = [1 if j in x else 0 for x in rang]
    pp = s.index(1)
    ax.scatter(bp_c[i],rsid_c[i], color=pcol[pp],s=19,alpha=0.7,linewidths=3)

ax.set_title(r'Manhattan plot for genetic variants tested for trans effects ', fontsize=40,fontweight="bold")
ax.set_xlim(0,2900)
ax.set_ylim(0, 11)
#ax.set_xlabel("Chromosome wide",{'size': axis_font_size, 'color': bordercolor},labelpad = 25)
ax.set_ylabel(r"$-\log_{10}(Pval)$",{'size': axis_font_size+10, 'color': bordercolor},labelpad = 25, fontweight="bold")
ax.plot(np.arange(0,2910,30),np.repeat(7.3,np.arange(0,2910,30).shape[0]),'r--',color="black",lw=2)
ax.plot(np.arange(0,2910,30),np.repeat(4.3,np.arange(0,2910,30).shape[0]),'r--',color="black",lw=2)
for side, border in ax.spines.items():
    border.set_linewidth(borderwidth)
    border.set_color(bordercolor)
plt.savefig("Manhattanplot_TranseQTLs.png")
plt.clf()














