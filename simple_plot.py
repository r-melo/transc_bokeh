from bokeh.io import show, output_notebook, curdoc
from bokeh.layouts import row, column, widgetbox, gridplot
from bokeh.models import ColumnDataSource, RangeTool, CustomJS, Select
from bokeh.plotting import figure, curdoc
from bokeh.models import  Band, HoverTool, Legend, Spacer, Slider, RangeSlider, ToolbarBase
from bokeh.palettes import Category10, Category20b, Category20c
from bokeh.models.widgets import CheckboxGroup, Panel, Tabs, Button
from bokeh.models.widgets import DataTable, DateFormatter, TableColumn
from bokeh.server.server import Server

from os.path import dirname, join
from difflib import SequenceMatcher
import numpy as np
import pandas as pd
import csv
import sys
import copy


# PVAL OUT OF RANGE ??



global argv
global lines
global pvals
global bands
global go_lines
global go_bands
global gocheckbox_group

def update(atrr, old, new):
    for i in range(len(lines)):
        if i in new:
            lines[i].visible = True
            bands[i].visible = True
            if (len(lines)) not in new:
                 bands[i].visible = False
        else:
            lines[i].visible = False
            bands[i].visible = False

        if (len(lines)) not in new:
            for i in range(len(lines)):
                bands[i].visible = False
                if i in new:
                    lines[i].visible = True

def update_transc(atrr, old, new):
    for i in range(len(lines)):
        if i in new:
            lines[i].visible = True
        else:
            lines[i].visible = False

def update_pval(atrr, old, new):
    for i in range(len(lines)):
        if i in new:
            lines[i].visible = True
            pvals[i].visible = True
            bands[i].visible = True
            if (len(lines)) not in new:
                 bands[i].visible = False
        else:
            lines[i].visible = False
            bands[i].visible = False
            pvals[i].visible = False

        if (len(lines)) not in new:
            for i in range(len(lines)):
                bands[i].visible = False
                if i in new:
                    lines[i].visible = True

                    
def update_go(atrr, old, new):
    for i in range(len(go_lines)):
        if i in new:
            go_lines[i].visible = True
            go_bands[i].visible = True
        else:
            go_lines[i].visible = False
            go_bands[i].visible = False

def up_sel():
    gocheckbox_group.active = [i for i in range(len(go_lines))]
    for i in range(len(go_lines)):
        go_lines[i].visible = True
        go_bands[i].visible = True
        
def up_uns():
    gocheckbox_group.active = []
    for i in range(len(go_lines)):
        go_lines[i].visible = False
        go_bands[i].visible = False
        
        
def plot(doc):
             
    global argv
    global lines
    global pvals
    global bands
    global go_lines
    global go_bands
    global gocheckbox_group

#    argl = len(argv)
    signal = argv[1]
    av_name = argv[2] 
    if argv[3] != 'NONE':
        pv_name = argv[3] 
    if argv[4] != 'NONE':
        go_name = argv[4]
    if argv[5] != 'NONE':
        kg_name = argv[5] 
    if argv[6] != 'NONE':
        dc_name = argv[6]

    lineskip=None
    transfile=False
    with open(av_name) as f:
        line = f.readline()
        if line[0] == '#':
            lineskip = 1
            transfile=True
            
    avreldf = pd.read_csv(av_name, sep='\t',skiprows=lineskip)
    
    ymax = 0
    ymin = 100000000000
    for col in avreldf.columns[2:]:
        if avreldf[col].max() > ymax:
            ymax = avreldf[col].max()
        if avreldf[col].min() < ymin:
            ymin = avreldf[col].min()
    ymin = min(0,ymin)
    
    colnames = ['Protein', 'Abs. Position']
    if not transfile:
        catsize = []
        catnames = []
        itercab = iter(avreldf.columns[2:])
        for indcat in range(len(avreldf.columns[2:])//2):
            cat = next(itercab)
            catname = cat.split()
            catsize += [len(catname) -2]
            if len(catname) > 3:
                string1 = catname[1]
                string2 = catname[2]
                match = SequenceMatcher(None, string1, string2).find_longest_match(1, len(string1)-1, 1, len(string2)-1)
                catnames += [string1[match.a: match.a + match.size]]
            else:
                catnames += catname[1][1:-1]
            cat = next(itercab)
        
        for i, col in enumerate(avreldf.iloc[:,3::2]):
            avreldf[col] = avreldf[col]/np.sqrt(catsize[i])
        
        for cat in catnames:
            colnames += ['Av. ' + cat]
            colnames += [cat]

        avreldf.columns = colnames
        
        avreldf['Position'] = avreldf['Abs. Position'] / avreldf['Abs. Position'].max()
        avreldf.sort_values('Position', inplace=True)
        for i, cat in enumerate(catnames):
            avreldf[cat + '_low'] = avreldf['Av. ' + cat] - avreldf[cat]
            avreldf[cat + '_up'] = avreldf['Av. ' + cat] + avreldf[cat]
            
    else:
        catnames = avreldf.columns[2:]
        for cat in catnames:
            colnames += ['Av. ' + cat]
        avreldf.columns = colnames
        avreldf['Position'] = avreldf['Abs. Position'] / avreldf['Abs. Position'].max()
        avreldf.sort_values('Position', inplace=True)

       
        
    
    if argv[3] != 'NONE':
        pvdf = pd.read_csv(pv_name, sep='\t')
        pvdf['Position'] = pvdf['position(dim1)'] / pvdf['position(dim1)'].max()
        pvdf.sort_values('Position', inplace=True)
    
        avreldf['pvalue ' + catnames[0]] = 1
        for i, cat in enumerate(catnames[1:]):
            avreldf['pvalue ' + cat] = pvdf.iloc[:,2+i]
        
    
    n = len(catnames)
    
    n_old = len(avreldf.columns)
        
    if argv[4] != 'NONE':
        godf = pd.read_csv(go_name, sep='\t')
        godf.rename(columns={'Prot' : 'Protein', 'Pos': 'Abs. Position'}, inplace=True)
        
        
        avreldf = pd.merge(avreldf,godf, how='left', left_on=['Protein', 'Abs. Position'], right_on=['Protein', 'Abs. Position'])
        
    if argv[5] != 'NONE':
        kgdf = pd.read_csv(kg_name, sep='\t')
        kgdf.rename(columns={'Prot' : 'Protein', 'Pos': 'Abs. Position'}, inplace=True)
            
        avreldf = pd.merge(avreldf,kgdf, how='left', left_on=['Protein', 'Abs. Position'], right_on=['Protein', 'Abs. Position'])
        
    n_new = len(avreldf.columns)
        
        
            
    avreldf['zero'] = 0
    
    

    tools_to_show = 'box_zoom,pan,save,hover,reset'
    ar = figure(plot_height =300, plot_width = 1200, 
                x_range=(0,1),
                y_axis_label = signal,
                tools=tools_to_show,
                toolbar_location='above')
    
    if argv[3] != 'NONE' or argv[3] != 'NONE':
        ar.xaxis.major_tick_line_color = None 
        ar.xaxis.minor_tick_line_color = None 
        ar.xaxis.major_label_text_font_size = '0pt'
    
    if argv[3] != 'NONE':
        pv = figure(plot_height =300, plot_width = 1200, 
                    toolbar_location=None,
                    tools='pan',
                    y_axis_label = 'p_values',
                    x_range=ar.x_range,
                    y_axis_type="log",
                    y_range=(10**(-5),1))    
        
    if argv[4] != 'NONE' or argv[5] != 'NONE':
        pv.xaxis.major_tick_line_color = None 
        pv.xaxis.minor_tick_line_color = None 
        pv.xaxis.major_label_text_font_size = '0pt'  
        
        go = figure(plot_height =300, plot_width = 1200, 
                    toolbar_location=None,
                    tools='hover,box_zoom',
                    y_axis_label = 'Ontologies',
                    x_range=ar.x_range)
        
    c_palette = ['#000000'] + Category10[10]
    legend_it = []
    glyph_list = []
    

    source = ColumnDataSource(avreldf)
    
    
    lines = []
    bands = []
    
    for i, cat in enumerate(catnames):
        if not transfile:
            bands.append(Band(base='Position', lower=cat + '_low', upper=cat + '_up', source=source, level='underlay', fill_alpha=0.3, line_alpha=0, fill_color=c_palette[i]))
            ar.add_layout(bands[i])
            
        lines.append(ar.line('Position', 'Av. ' + cat, source=source, color=c_palette[i % 11], name=cat))
        
        legend_it.append((cat, [lines[i]]))
        
    legend = Legend(items=legend_it)
    
                 
    ar.add_layout(legend, 'right')
    
    hover = ar.select(dict(type=HoverTool))
    if transfile:
        hover.tooltips = [ ("Group", "$name"),("Position", "@{Abs. Position}"), ("Protein", "@Protein")]
    else:
        hover.tooltips = [ ("Group", "$name"),("Position", "@{Abs. Position}"), ("Protein", "@Protein"),
                          ("Std.", "@$name")]
    hover.mode = 'mouse'
    
    if argv[3] != 'NONE':
        pvals = []
        
        for i, cat in enumerate(catnames):
            pvals.append(pv.line('Position', 'pvalue ' + cat, source=source, color = c_palette[i], line_width=0.5))
        
    
    if argv[4] != 'NONE' or argv[5] != 'NONE':
        go_palette = Category20b[20] + Category20c[20]
        go_bands = []
        go_lines = []
        go_legends = []
        gonames = []
        for i, col in enumerate(avreldf.columns[n_old:n_new]): #:n_new].columns):
            gonames.append(col)
            go_bands.append(Band(base='Position', upper=col, lower='zero', source=source, fill_alpha=0.5, line_alpha=0.2, fill_color=go_palette[i % 40])) #, muted_color=go_palette[i % 40],  muted_alpha=0.2, muted=True
            go.add_layout(go_bands[i])
            go_lines.append(go.line('Position', col, source=source, line_alpha=0.7, color=go_palette[i % 40], name=col))
            #print(i,col)
            go_legends.append((col, [go_lines[i]]))
                 
        golegend = Legend(items=go_legends)

        gohover = go.select(dict(type=HoverTool))
        gohover.tooltips = [ ("Ontology", "$name")]
        hover.mode = 'mouse'
        
    
    select = []
    callbacks = []
    cnames = ['Black', 'Blue', 'Orange', 'Green', 'Red', 'Purple', 'Brown', 'Pink', 'Grey', 'Olive', 'Cyan']
    cdict = list(zip(c_palette,cnames))
    
    print(argv[3])
    print(argv[4])
    if argv[3] == 'NONE':
        if transfile:
            for i in range(len(lines)):
                select.append(Select(title=catnames[i] + " color", value=' ', options=cdict))
                callbacks.append(CustomJS(args=dict(source=source, line=lines[i], select=select[i], dic=cdict), code ="""
                line.glyph.line_color = select.value;
                line.trigger('change');
                """))
                select[i].js_on_change('value', callbacks[i])
        else:
            for i in range(len(lines)):
                select.append(Select(title=catnames[i] + " color", value=' ', options=cdict))
                callbacks.append(CustomJS(args=dict(source=source, line=lines[i], band=bands[i], select=select[i], dic=cdict), code ="""
                line.glyph.line_color = select.value;
                band.fill_color = select.value;
                line.trigger('change');
                band.trigger('change');
                """))
            
                select[i].js_on_change('value', callbacks[i])
    
    else:
        for i in range(len(lines)):
            select.append(Select(title=catnames[i] + " color", value=' ', options=cdict))
            callbacks.append(CustomJS(args=dict(source=source, line=lines[i], band=bands[i], pval=pvals[i], select=select[i], dic=cdict), code ="""
            pval.glyph.line_color = select.value;
            line.glyph.line_color = select.value;
            band.fill_color = select.value;
            line.trigger('change');
            band.trigger('change');
            pval.trigger('change')
            """))
    
            select[i].js_on_change('value', callbacks[i])
    
    if argv[4] == 'NONE' and argv[5] == 'NONE':
        if transfile:
            checkbox_group = CheckboxGroup(labels = list(catnames), active = [i for i in range(len(catnames)+1)])
            checkbox_group.on_change('active',update_transc)
        else:
            checkbox_group = CheckboxGroup(labels = list(catnames + ["Std."]), active = [i for i in range(len(catnames)+1)])
            checkbox_group.on_change('active',update)
    else:
        checkbox_group = CheckboxGroup(labels = list(catnames + ["Std."]), active = [i for i in range(len(catnames)+1)])
        checkbox_group.on_change('active',update_pval)
    
    if argv[4] != 'NONE' or argv[5] != 'NONE':
        gocheckbox_group = CheckboxGroup(labels = list(gonames), active = [i for i in range(len(gonames))])
        gocheckbox_group.on_change('active',update_go)
        
        sel_all = Button(label="Select All")
        uns_all = Button(label="Unselect All")
        
        sel_all.on_click(up_sel)
        uns_all.on_click(up_uns)
    
    if argv[6] == 'NONE':
        tabdf = avreldf[['Protein', 'Abs. Position']]
        maxim=tabdf['Abs. Position'].max()
        tabdf = tabdf.set_index(tabdf['Abs. Position'].values)
        tabdf = tabdf.reindex([i for i in range(maxim)], fill_value=' ')
        tabdf.to_csv('toto.txt', sep='\t')
        tab_source = ColumnDataSource({'Protein' : [], 'Abs. Position' : []})
        orig_source = ColumnDataSource({'Protein' : [], 'Abs. Position' : []})
        tab_source.data = {
            'protein' : tabdf.Protein,
            'position' : tabdf['Abs. Position']}
        orig_source.data = {
            'protein' : tabdf.Protein,
            'position' : tabdf['Abs. Position']}
        columns = [
                TableColumn(field="protein", title="Protein"),
                TableColumn(field="position", title="Position"),
            ]
        data_table = DataTable(source=tab_source, columns=columns, width=300, height=600)
    
        button = Button(label="Download", button_type="success")
        p_button = Button(label="Proteins Only", button_type="success")
    
        button.callback = CustomJS(args=dict(source=tab_source),code="""
            var data = source.data;
            var filetext = 'protein,position\\n';
            for (var i = 0; i < data['protein'].length; i++) {
                var currRow = [data['protein'][i].toString(),
                               data['position'][i].toString().concat('\\n')];
            
                var joined = currRow.join();
                filetext = filetext.concat(joined);
            }
            
            var filename = 'data_result.csv';
            var blob = new Blob([filetext], { type: 'text/csv;charset=utf-8;' });
            
            //addresses IE
            if (navigator.msSaveBlob) {
                navigator.msSaveBlob(blob, filename);
            } else {
                var link = document.createElement("a");
                link = document.createElement('a')
                link.href = URL.createObjectURL(blob);
                link.download = filename
                link.target = "_blank";
                link.style.visibility = 'hidden';
                link.dispatchEvent(new MouseEvent('click'))
            }
            """)

        p_button.callback = CustomJS(args=dict(source=tab_source),code="""
            var data = source.data;
            var filetext = 'protein\\n';
            for (var i = 0; i < data['protein'].length; i++) {
                var currRow = [data['protein'][i].toString().concat('\\n')];
            
                var joined = currRow.join();
                filetext = filetext.concat(joined);
            }
            
            var filename = 'data_result.csv';
            var blob = new Blob([filetext], { type: 'text/csv;charset=utf-8;' });
            
            //addresses IE
            if (navigator.msSaveBlob) {
                navigator.msSaveBlob(blob, filename);
            } else {
                var link = document.createElement("a");
                link = document.createElement('a')
                link.href = URL.createObjectURL(blob);
                link.download = filename
                link.target = "_blank";
                link.style.visibility = 'hidden';
                link.dispatchEvent(new MouseEvent('click'))
            }
            """)

    
    else:
        dictdf = pd.read_csv(dc_name, sep='\t',skiprows=lineskip, names=['A','B'])
        dictdf['A'] = dictdf['A'].astype(str)
        dictdf['B'] = dictdf['B'].astype(str)
        dictdf.drop_duplicates(subset='A', inplace=True)
        dictdf.drop_duplicates(subset='B', inplace=True)
        tabdf = avreldf[['Protein', 'Abs. Position']]
        
        if not set(tabdf['Protein'].values).isdisjoint(set(map(lambda x:x.upper(),dictdf['A'].values))):
            dictdf['A'] = dictdf['A'].str.upper()
            tabdf = pd.merge(tabdf,dictdf, how='left', left_on=['Protein'], right_on=['A'])
            tabdf.rename(columns={'B' : 'Translated'}, inplace=True)
        if not set(tabdf['Protein'].values).isdisjoint(set(map(lambda x:x.lower(),dictdf['A'].values))):
            dictdf['A'] = dictdf['A'].str.lower()
            tabdf = pd.merge(tabdf,dictdf, how='left', left_on=['Protein'], right_on=['A'])
            tabdf.rename(columns={'B' : 'Translated'}, inplace=True)
        if not set(tabdf['Protein'].values).isdisjoint(set(map(lambda x:x.upper(),dictdf['B'].values))):
            dictdf['B'] = dictdf['B'].str.upper()
            tabdf = pd.merge(tabdf,dictdf, how='left', left_on=['Protein'], right_on=['B'])
            tabdf.rename(columns={'A' : 'Translated'}, inplace=True)
        if not set(tabdf['Protein'].values).isdisjoint(set(map(lambda x:x.lower(),dictdf['B'].values))):
            dictdf['B'] = dictdf['B'].str.lower()
            tabdf = pd.merge(tabdf,dictdf, how='left', left_on=['Protein'], right_on=['B'])
            tabdf.rename(columns={'A' : 'Translated'}, inplace=True)
        
        
        maxim=tabdf['Abs. Position'].max()
        tabdf = tabdf.set_index(tabdf['Abs. Position'].values)
        tabdf = tabdf.reindex([i for i in range(maxim)], fill_value=' ')
        tabdf.to_csv('toto.txt', sep='\t')
        tab_source = ColumnDataSource({'Protein' : [], 'Abs. Position' : [], 'Translated' : []})
        orig_source = ColumnDataSource({'Protein' : [], 'Abs. Position' : [], 'Translated' : []})
        tab_source.data = {
            'protein' : tabdf.Protein,
            'position' : tabdf['Abs. Position'],
            'translated' : tabdf.Translated}  
        orig_source.data = {
            'protein' : tabdf.Protein,
            'position' : tabdf['Abs. Position'],
            'translated' : tabdf.Translated}
        columns = [
                TableColumn(field="protein", title="Protein"),
                TableColumn(field="position", title="Position"),
                TableColumn(field="translated", title="Translated"),
            ]
        data_table = DataTable(source=tab_source, columns=columns, width=300, height=600)
    
        button = Button(label="Download", button_type="success")
        p_button = Button(label="Proteins Only", button_type="success")
        t_button = Button(label="Translated Only", button_type="success")
        
        button.callback = CustomJS(args=dict(source=tab_source),code="""
            var data = source.data;
            var filetext = 'protein,position,translated\\n';
            for (var i = 0; i < data['protein'].length; i++) {
                var currRow = [data['protein'][i].toString(),
                               data['position'][i].toString(),
                               data['translated'][i].concat('\\n')];
            
                var joined = currRow.join();
                filetext = filetext.concat(joined);
            }
            
            var filename = 'data_result.csv';
            var blob = new Blob([filetext], { type: 'text/csv;charset=utf-8;' });
            
            //addresses IE
            if (navigator.msSaveBlob) {
                navigator.msSaveBlob(blob, filename);
            } else {
                var link = document.createElement("a");
                link = document.createElement('a')
                link.href = URL.createObjectURL(blob);
                link.download = filename
                link.target = "_blank";
                link.style.visibility = 'hidden';
                link.dispatchEvent(new MouseEvent('click'))
            }
            """)

        p_button.callback = CustomJS(args=dict(source=tab_source),code="""
            var data = source.data;
            var filetext = 'protein\\n';
            for (var i = 0; i < data['protein'].length; i++) {
                var currRow = [data['protein'][i].toString().concat('\\n')];
            
                var joined = currRow.join();
                filetext = filetext.concat(joined);
            }
            
            var filename = 'data_result.csv';
            var blob = new Blob([filetext], { type: 'text/csv;charset=utf-8;' });
            
            //addresses IE
            if (navigator.msSaveBlob) {
                navigator.msSaveBlob(blob, filename);
            } else {
                var link = document.createElement("a");
                link = document.createElement('a')
                link.href = URL.createObjectURL(blob);
                link.download = filename
                link.target = "_blank";
                link.style.visibility = 'hidden';
                link.dispatchEvent(new MouseEvent('click'))
            }
            """)
        
        t_button.callback = CustomJS(args=dict(source=tab_source),code="""
            var data = source.data;
            var filetext = 'protein\\n';
            for (var i = 0; i < data['translated'].length; i++) {
                var currRow = [data['translated'][i].toString().concat('\\n')];
            
                var joined = currRow.join();
                filetext = filetext.concat(joined);
            }
            
            var filename = 'data_result.csv';
            var blob = new Blob([filetext], { type: 'text/csv;charset=utf-8;' });
            
            //addresses IE
            if (navigator.msSaveBlob) {
                navigator.msSaveBlob(blob, filename);
            } else {
                var link = document.createElement("a");
                link = document.createElement('a')
                link.href = URL.createObjectURL(blob);
                link.download = filename
                link.target = "_blank";
                link.style.visibility = 'hidden';
                link.dispatchEvent(new MouseEvent('click'))
            }
            """)
        
        
        
    
    y_slider_callback = CustomJS(args=dict(y_range=ar.y_range),code="""
        var range = cb_obj.value;
        y_range.setv({"start": range[0], "end": range[1]});
    """)
    y_slider = RangeSlider(title=signal + " Axis Range", start=ymin, end=int(-(-ymax // 1)), value=(ymin+0.2, ymax-0.2), step=.2, callback=y_slider_callback)
    y_slider.js_on_change('value', y_slider_callback)
    
    
    
    ar.x_range.callback = CustomJS(args=dict(source=tab_source, orig=orig_source, maxim=maxim),code="""
        var start = cb_obj.start;
        var end = cb_obj.end;
        var d1 = orig.data;
        var d2 = source.data;
        start *= maxim;
        end *= maxim;
        start = parseInt(start) + 1;
        end = parseInt(end) + 1;
        d2['protein'] = d1['protein'].slice(start,end);
        d2['position'] = d1['position'].slice(start,end);
        source.change.emit();
    """)

    
        
    
        

    if argv[4] != 'NONE' or argv[5] != 'NONE':
        ar_tab = Panel(child=gridplot([[ar],[pv],[go]], toolbar_location='above', merge_tools=True), title=signal)
        go_tab = Panel(child=row(gocheckbox_group,widgetbox(sel_all,uns_all)), title='Ontologies')
    elif argv[3] != 'NONE':
        ar_tab = Panel(child=gridplot([[ar],[pv]], toolbar_location='above', merge_tools=True), title=signal)
    else:
        ar_tab = Panel(child=ar, title=signal)
    
    ctrl_tab = Panel(child=row(checkbox_group,widgetbox(select),widgetbox(y_slider)), title='Controls')
    if argv[6] != 'NONE':
        table_tab = Panel(child=row(widgetbox(data_table),widgetbox(button,p_button,t_button)), title='Proteins')
    else:
        table_tab = Panel(child=row(widgetbox(data_table),widgetbox(button,p_button)), title='Proteins')
    if argv[4] != 'NONE' or argv[5] != 'NONE':
        tabs = Tabs(tabs=[ar_tab, ctrl_tab, table_tab, go_tab])
    else:
        tabs = Tabs(tabs=[ar_tab, ctrl_tab, table_tab])
    #tabs.show()
    doc.add_root(tabs)

argv = sys.argv
    
server = Server({'/': plot}, num_procs=1)
server.start()

server.io_loop.add_callback(server.show, "/")
server.io_loop.start()

