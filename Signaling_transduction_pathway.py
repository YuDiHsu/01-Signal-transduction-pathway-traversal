import datetime
import os
import pickle
import re
from collections import deque
from multiprocessing import freeze_support

import pandas as pd
import requests
import xlsxwriter
from bs4 import BeautifulSoup
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import WebDriverWait


def _balanced_bracket(t):
    count = 0
    # print(t)
    n = len(t)
    new_t = ''
    # Maintain a count for opening
    # brackets Traversing string
    for i in range(n):
        # # check if opening bracket
        if t[i] == '(':
            # # print str[i] and increment
            # # count by 1
            new_t += t[i]
            # # print(t[i], end="")
            count += 1
        # # check if closing bracket and count != 0
        elif t[i] == ')' and count != 0:
            new_t += t[i]
            # # print(t[i], end="")
            # # decrement count by 1
            count -= 1
        # # if str[i] not a closing brackets
        # # print it
        elif t[i] != ')':
            new_t += t[i]
            # print(t[i], end="")
    # # balanced brackets if opening brackets
    # # are more then closing brackets
    if count != 0:
        # # print remaining closing brackets
        for j in range(count):
            new_t += ')'
            # print(")", end="")
    return new_t


def _parse_nested(text, left=r'[(]', right=r'[)]', sep=r','):
    # # check and balanced bracket
    bt = _balanced_bracket(text)
    pat = r'({}|{}|{})'.format(left, right, sep)
    tokens = re.split(pat, bt)
    stack = [[]]
    for x in tokens:
        if not x or re.match(sep, x):
            continue
        if re.match(left, x):
            stack[-1].append([])
            stack.append(stack[-1][-1])
        elif re.match(right, x):
            stack.pop()
            if not stack:
                raise ValueError('error: opening bracket is missing')
        else:
            _x = x.strip('+').rstrip('+')
            if _x:
                for _ in _x.split('+'):
                    stack[-1].append(_)
    if len(stack) > 1:
        print(stack)
        raise ValueError('error: closing bracket is missing')

    return stack.pop()[0]


def _route_processing(r_list):
    new_pathway_id_route_list = []
    # # inhibition position ex: A -> B -| C pos:2
    inhibit_idx_list = []
    # # involvement position ex: A -- B -| C pos:1
    involve_idx_list = []
    record_route_idx = 0
    for _route in r_list:
        if re.search(' -\| | =\| ', _route):
            sub_route = re.split(' -\| | =\| ', _route)
            for sr in sub_route:
                if re.search(' -- ', sr):
                    sub_route_2 = re.split(' -- ', sr)
                    new_pathway_id_route_list += [_parse_nested(i) for i in sub_route_2]
                    involve_idx_list.append(len(new_pathway_id_route_list) - 1)
                    record_route_idx = len(new_pathway_id_route_list) - 1
                else:
                    new_pathway_id_route_list.append(_parse_nested(sr))
                    record_route_idx = len(new_pathway_id_route_list) - 1
            if record_route_idx:
                inhibit_idx_list.append(record_route_idx)

        else:
            if re.search(' -- ', _route):
                sub_route_2 = re.split(' -- ', _route)
                new_pathway_id_route_list += [_parse_nested(_i) for _i in sub_route_2]
                if len(new_pathway_id_route_list) - 1 not in involve_idx_list:
                    involve_idx_list.append(len(new_pathway_id_route_list) - 1)
            else:
                # print('r_route', _route)
                new_pathway_id_route_list.append(_parse_nested(_route))

    print(new_pathway_id_route_list)
    # if inhibit_idx_list:
    #     print('inhibit at: ', inhibit_idx_list)
    # if involve_idx_list:
    #     print('involve at', involve_idx_list)
    return [new_pathway_id_route_list, inhibit_idx_list, involve_idx_list]


def _cross_combine(arr_1):
    if isinstance(arr_1, list):
        d_1 = deque(arr_1)
    else:
        d_1 = [arr_1]
    result_1 = []
    while len(d_1):
        e = d_1.pop()
        if isinstance(e, str):
            result_1.append(e)
        elif isinstance(e, list):
            for _e in e:
                d_1.append(_e)

    return result_1


def _involved_sym_id(data_dict_array: dict):
    new_data_dict_array = []
    for data in data_dict_array:
        gene_id = data['GeneID'].replace('hsa:', '').split('/')
        data['GeneIDFlatten'] = gene_id
        symbol = data['Symbol'].split('/')
        data['SymbolFlatten'] = symbol
        data['IDpairSym'] = {}
        try:
            for idx, _id in enumerate(gene_id):
                if _id not in data['IDpairSym']:
                    data['IDpairSym'][_id] = symbol[idx]
                else:
                    print('exception: ', symbol)
                    print('exception: ', _id)
                    print('temp: ', data['IDpairSym'][_id])

        except ValueError as V:
            print(V)
            print(gene_id, '\t', (len(gene_id)))
            print(symbol, '\t', (len(symbol)))
        new_data_dict_array.append(data)

    return new_data_dict_array


def _selenium_parser_network_id():
    options = webdriver.ChromeOptions()
    options.add_argument('--headless')
    options.add_argument('--disable-gpu')

    # 初始化selenium driver時傳入option參數
    driver = webdriver.Chrome(chrome_options=options)
    driver.get("https://www.kegg.jp/kegg-bin/search_pathway_text?map=hsa&keyword=&mode=1&viewImage=false")

    network_id_list = []
    hsaMap_network_dict = {}
    for i in range(18):
        WebDriverWait(driver, 5).until(
            EC.presence_of_element_located((By.XPATH, '/html/body/form/input[12]')))
        if i == 0:
            driver.find_element_by_xpath('/html/body/form/input[12]')
        else:
            driver.find_element_by_xpath('/html/body/form/input[12]').click()

        soup = BeautifulSoup(driver.page_source, "html.parser")
        html_content = soup.find_all("a")
        for link in html_content:
            if link.text not in hsaMap_network_dict:
                hsaMap_network_dict[link.text] = []
            if re.search('hsa\d+', link.text):
                _url = f'https://www.genome.jp/pathway/{link.text}'
                _res = requests.get(_url)
                _soup = BeautifulSoup(_res.text, 'html.parser')
                cont = _soup.findAll('li', {'class': 'entry class-D network'})
                print('--', link.text, '--')
                for t in cont:
                    nid = t.text.strip().split('\n\t\t\n\t\t')[0]
                    hsaMap_network_dict[link.text].append(nid)
                    if nid not in network_id_list:
                        network_id_list.append(nid)
                        print(nid)

    with open(os.path.join('.', 'pickle_storage', 'total_networkID_list.pickle'), 'wb') as r:
        pickle.dump(network_id_list, r)
    with open(os.path.join('.', 'pickle_storage', 'total_hsaMap_network_dict.pickle'), 'wb') as r:
        pickle.dump(hsaMap_network_dict, r)

    # return os.path.join('.', 'pickle_storage', 'total_networkID_list.pickle')


def _bs4_parser_info_from_kegg():
    with open(os.path.join('.', 'pickle_storage', 'total_networkID_list.pickle'), 'rb') as r:
        network_id_list = pickle.load(r)

    run_times = 0
    pathway_info_list = []
    counter = 0

    # # analyze and transform info to acceptable format
    try:
        while run_times < 3:
            for i in network_id_list:
                sub_pathway_url = f'https://www.genome.jp/entry/{i}'
                temp_df = pd.read_html(sub_pathway_url)

                # # check in hsa pathway
                check_pathway_col = temp_df[3][0].values.tolist()
                if 'Pathway' in check_pathway_col:
                    check_pathway_name = temp_df[3][1][temp_df[3][0].values.tolist().index('Pathway')]
                    if re.search('hsa\d+', check_pathway_name):
                        print(sub_pathway_url)
                        subPathwayID = temp_df[3][1][0].replace(' Network', '')
                        subPathway = temp_df[3][1][1]
                        raw_symbol_route = temp_df[3][1][2]
                        raw_id_route = temp_df[3][1][3]

                        # # rearrange order of pathway_id_route_list
                        symbol_route_list = _route_processing(re.split(' -> | => | >> | // | -- ', raw_symbol_route))
                        id_route_list = _route_processing(re.split(' -> | => | >> | // | -- ', raw_id_route))

                        # # symbol and geneID flatten to List
                        symbol_flatten_list = _cross_combine(symbol_route_list[0])
                        geneID_flatten_list = _cross_combine(id_route_list[0])

                        # m_p = temp_df[3][1][temp_df[3][0].values.tolist().index('Pathway')]
                        re_m_p = [ele.strip() for ele in re.split('hsa\d+\s', check_pathway_name) if ele]
                        id_m_p = re.findall('hsa\d+', check_pathway_name)
                        try:
                            for id_ in id_m_p:
                                # # Store info to dict
                                template_dict = {'metaPathway': re_m_p[id_m_p.index(id_)], 'metaID': id_,
                                                 'subPathway': subPathway,
                                                 'subPathwayID': subPathwayID, 'symRoute': raw_symbol_route,
                                                 'involvedGeneID': geneID_flatten_list[::-1],
                                                 'involvedSymID': symbol_flatten_list[::-1]}

                                if template_dict not in pathway_info_list:
                                    pathway_info_list.append(template_dict)
                                    counter += 1
                                    print(template_dict)
                                    print(f'Data Processing successfully complete, length = {counter}')
                                    print('=' * 100)
                        except IndexError as Idx:
                            if len(re_m_p) != len(id_m_p):
                                print(print(sub_pathway_url))
                                print(f'Length is different, name_length: {re_m_p}, id_length: {id_m_p}')
                                print('*' * 100)
            break

    except RuntimeError as R:
        run_times += 1
        print(R)

    # # dump data to pickle
    with open(os.path.join('.', 'pickle_storage', 'total_pathway_info_list.pickle'), 'wb') as f:
        pickle.dump(pathway_info_list, f)

    # # # validate data
    # with open(os.path.join('.', 'pickle_storage', 'total_pathway_info_list.pickle'), 'rb') as r:
    #     networkInfo_list = pickle.load(r)


def _data_processing(new_data, all_info_data):
    analyzed_data_list = []
    for n_d in new_data:
        for info in all_info_data:
            involved_sym_name = []
            total_len = 0
            flag = False
            re_n_d = {}
            for id_ in n_d['GeneIDFlatten']:
                if id_ in info['involvedGeneID']:
                    total_len = len(info['involvedGeneID'])
                    flag = True
                    if n_d['IDpairSym'][id_] not in involved_sym_name:
                        involved_sym_name.append(n_d['IDpairSym'][id_])
            if flag:
                re_n_d['metaPathway'] = info['metaPathway']
                re_n_d['metaID'] = info['metaID']
                re_n_d['subPathway'] = info['subPathway']

                if total_len and involved_sym_name:
                    re_n_d['Ratio'] = round(len(involved_sym_name) / total_len, 2)
                else:
                    re_n_d['Ratio'] = 0
                if involved_sym_name:
                    re_n_d['InvolvedSym'] = '/'.join(involved_sym_name)
                else:
                    re_n_d['InvolvedSym'] = ''

                re_n_d['SymRoute'] = info['symRoute']
                re_n_d['subPathwayID'] = info['subPathwayID']
                re_n_d['GeneID'] = n_d['GeneID']
                re_n_d['RNASeqSym'] = n_d['Symbol']
                if re_n_d not in analyzed_data_list:
                    analyzed_data_list.append(re_n_d)
    return analyzed_data_list


def _write_to_xlsx(sig_data: list, sh_name: list, header_name_list: list):
    date = datetime.datetime.now().date()
    workbook = xlsxwriter.Workbook(
        os.path.join('.', 'exported_data', f"PathwayTraversal_{date.strftime('%Y%m%d')}.xlsx"))
    col_len_width = []
    for idx, data_ in enumerate(zip(sig_data, sh_name)):
        header_format = workbook.add_format(
            {'align': 'center', 'valign': 'vcenter', 'size': 12, 'color': 'black', 'bold': 4})
        header = []
        for h in header_name_list:
            header.append({'header': h, "format": header_format})

        # print(sheet_name)
        sheet = workbook.add_worksheet(data_[1])
        if len(data_[0]):

            sheet.add_table(0, 0, len(data_[0]), len(data_[0][0]) - 1,
                            {'data': sorted(data_[0], key=lambda x: x[3], reverse=True), 'autofilter': True,
                             'columns': header})

            for m in range(len(data_[0]) + 1):
                sheet.set_row(m, 30, cell_format=header_format)

            col_len_width.append([len(j) for j in header_name_list])
            for n, l in enumerate(zip(col_len_width[idx], header_name_list)):
                sheet.set_column(n, n, max(l[0], len(l[1])) * 5)
    workbook.close()

    return os.path.join('.', 'exported_data', f"PathwayTraversal_{date.strftime('%Y%m%d')}.xlsx")


def _traversal(raw_data_path):
    # # get raw data
    # # cols = ['ID', 'Description', 'GeneRatio.Item', 'GeneRatio.List', 'BgRatio.Item',
    #          'BgRatio.List', 'pvalue', 'p.adjust', 'qvalue', 'Count', 'GeneID',
    #          'Symbol', 'Datalink']
    df = pd.read_csv(raw_data_path, sep='\t')
    df = df.loc[:, ('ID', 'Description', 'GeneID', 'Symbol', 'Datalink')]
    data_dict_list = df.to_dict('records')

    # # involve data collection
    new_data_dict_list = _involved_sym_id(data_dict_list)

    # # load data to pickle
    with open(os.path.join('.', 'pickle_storage/total_pathway_info_list.pickle'), 'rb') as r:
        info_data_list = pickle.load(r)

    final_data = _data_processing(new_data_dict_list, info_data_list)

    # #　remove cols of RNASeqSym and GeneID to reduce duplicate
    final_df = pd.DataFrame(final_data)
    final_df = final_df.drop_duplicates(subset=['metaPathway', 'metaID', 'subPathway', 'Ratio', 'InvolvedSym',
                                                'SymRoute', 'subPathwayID'])
    # # keep cols of RNASeqSym and GeneID
    # final_df = pd.DataFrame(final_data)

    # # # multiprocessing
    # # # initiate multiprocessing
    # freeze_support()
    #
    # num_processes = int(mp.cpu_count() / 2)
    # chunk_size = int(len(new_data_dict_list) / num_processes)
    # chunks = [new_data_dict_list[i:i + chunk_size] for i in range(0, len(new_data_dict_list), chunk_size)]
    # pool = mp.Pool(processes=num_processes)
    # result_list = pool.map(partial(_data_processing, all_info_data=info_data_list), chunks)
    # final_data_list = []
    # for i in result_list:
    #     final_data_list += i
    # print(len(final_data_list))

    # # write to excel
    final_data_list = final_df.values.tolist()
    max_len = list(set([len(l) for l in final_df.values.tolist()]))
    _write_to_xlsx([final_data_list], ['Pathway_Traversal'], list(final_df.columns))


if __name__ == '__main__':
    # # make dir roots
    dir_roots = [os.path.join('.', 'imported_data'), os.path.join('.', 'pickle_storage')]
    for dr in dir_roots:
        if not os.path.exists(dr):
            os.makedirs(dr)

    # # parse network_id from kegg
    # _selenium_parser_network_id()

    # # parse info of pathway from kegg
    # _bs4_parser_info_from_kegg()

    # # traversal with RNAseq or Omics or merged data
    data_path = os.path.join('.', 'imported_data', 'MR.M_Ctrl_kegg_DEGAll.txt')
