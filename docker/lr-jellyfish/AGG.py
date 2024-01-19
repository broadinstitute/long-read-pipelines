import gzip
import json
import numpy

class GraphicalGenome:
    def __init__(self, filename):
        self.anchor, self.edges, self.outgoing, self.incoming = self.loadgraph(filename)

    def loadgraph(self, filename):
        if filename.endswith(".gz"):
            with gzip.open(filename, 'rb') as fp:
                data = fp.readlines()
        else:
            with open(filename, 'rb') as fp:
                data = fp.readlines()
        Anchor_dict = {}
        Edge_dict = {}
        Outgoing = {}
        Incoming = {}

        for line in data:
            line = line.decode()[:-1]
            if line.startswith("S"):
                itemlist = line.split('\t')
                name = itemlist[1]
                seq = itemlist[2]
                annotation = itemlist[3][5:]
                if name.startswith("A"):
                    Anchor_dict[name] = json.loads(annotation)
                    Anchor_dict[name]['seq'] = seq
                elif name.startswith("E"):
                    Edge_dict[name] = json.loads(annotation)
                    Edge_dict[name]['seq'] = seq
            elif line.startswith("L"):
                itemlist = line.split('\t')
                src = itemlist[1]
                dst = itemlist[3]
                Outgoing[src] = Outgoing.get(src, []) + [dst]
                Incoming[dst] = Incoming.get(dst, []) + [src]
                    
        return Anchor_dict, Edge_dict, Outgoing, Incoming
    
def reconstruct_refpath(graph, ref = 'NC_012920.1'):
    src = "SOURCE"
    Path = []

    visited = set()

    while src != "SINK":
        edgelist = graph.outgoing[src]
        ref_edge = []
        for edge in edgelist:
            if ref in graph.edges[edge]['reads']:
                ref_edge += [edge]
        if len(ref_edge) == 1:
            edge = ref_edge[0]
        else:
            for edge in ref_edge:
                if graph.outgoing[edge] == "SINK":
                    break
        Path += [src, edge]
        visited.add(src)
        src = graph.outgoing[edge][0]
        assert src.startswith("A") or src == "SINK"
        if src in visited:
            break

    reconstruct = ''
    for item in Path:
        if item == "SOURCE":
            reconstruct += ""
        elif item.startswith('A'):
            reconstruct += graph.anchor[item]['seq']
        elif item.startswith("E"):
            reconstruct += graph.edges[item]['seq']
        else:
            print(item)
    return reconstruct, Path

def write_gfa(AnchorInfo, Edge_info, outputfilename):
    header = ['H\tVN:Z:1.0\n']
    # add json annotation to edge and anchor sequences
    AnchorS = []
    for anchor,D in AnchorInfo.items():
        seq = D.pop('seq')
        json_string = json.dumps(D)
        AnchorS += ['S\t%s\t%s\t%s\n' % (anchor, seq, "PG:J:" + json_string)]

    EdgeS = []   
    Link = [] 
    for edge,edge_dict in Edge_info.items():
        seq = edge_dict.pop('seq')
        src = edge_dict.pop('src')
        dst = edge_dict.pop('dst')

        json_string = json.dumps(edge_dict)
        EdgeS += ['S\t%s\t%s\t%s\t%s\n' % (edge, seq, "PG:J:" + json_string, "RC:i:" + str(len(edge_dict['reads'])))]
        Link.append('L\t%s\t%s\t%s\t%s\t%s\n'% (src, "+", edge, "+", "0M"))
        Link.append('L\t%s\t%s\t%s\t%s\t%s\n'% (edge, "+", dst, "+", "0M"))     

    with open(outputfilename, 'w') as fp:
        for h in header:
            fp.write(h)
        for s in AnchorS:
            fp.write(s)
        for s in EdgeS:
            fp.write(s)
        for l in Link:
            fp.write(l)

def reconstruct_sample_subgraph(graph, samplelist, outputfile):
    Anchor_dict = {}
    Edge_dict = {}
    edgelist = list(graph.edges.keys())
    Anchor_dict = {}
    Edge_dict = {}
    nodelist = []
    for edge in edgelist:
        if len(set(samplelist) & set(graph.edges[edge]['strain'])) > 0:
            Edge_dict[edge] = graph.edges[edge]
            Edge_dict[edge]['src'] = graph.incoming[edge][0]
            Edge_dict[edge]['dst'] = graph.outgoing[edge][0]
        nodelist += graph.incoming[edge]
        nodelist += graph.outgoing[edge]

    nodelist = list(set(nodelist))
    for node in nodelist:
        if node.startswith("A"):
            Anchor_dict[node] = graph.anchor[node]

    write_gfa(Anchor_dict, Edge_dict, outputfile)

    return Anchor_dict, Edge_dict

def find_most_supported_path(graph, sample, ref = '012920'):
    src = sorted(graph.anchor.keys())[0]
    dst = sorted(graph.anchor.keys())[-1]
    Path = []

    visited = set()
    
    def find_most_supported_edge(graph, sample, src, ref):
        
        edgelist = graph.outgoing.get(src, [])
        sample_edge = []
        for edge in edgelist:
            if sample in graph.edges[edge]['strain']:
                sample_edge += [edge]
        if len(sample_edge)<1: # no haplotype for this sample at this place the go to the reference one
            for ref_edge in edgelist:
                if ref in set(graph.edges[ref_edge]['strain']):
                    break
            edge = ref_edge

        elif len(sample_edge) == 1:
            edge = sample_edge[0]
        else:
            # find the most supported edge
            read_count = [len(graph.edges[e]['reads']) for e in sample_edge]

            index = numpy.argmax(read_count)

            edge = sample_edge[index]
            
        return edge
            
    while src != dst:
        edge = find_most_supported_edge(graph, sample, src, ref)
        Path += [src, edge]

        visited.add(src)

        src = graph.outgoing[edge][0]

        assert src.startswith("A") or src == "SINK"

        if src in visited:
            break
    
    if src == dst:
        try:
            edge = find_most_supported_edge(graph, sample, src, ref)
            Path += [src, edge]
        except:
            Path = Path
        
    reconstruct = ''
    for item in Path:
        if item == "SOURCE":
            reconstruct += ""
        elif item.startswith('A'):
            reconstruct += graph.anchor[item]['seq']
        elif item.startswith("E"):
            reconstruct += graph.edges[item]['seq']
        else:
            print(item)

    return reconstruct, Path


class Find_all_Path_between_anchors:
    def __init__(self, graph, start, end, read_sets):
        self.subpath = []
        self.initial_set = read_sets
        self.find_path(graph, start, end, [], 0, self.initial_set)
        
    def find_path(self, g, start, end, sofar, depth, readset):
        
        if start == end:
            sofar1 = sofar + [end]
            if len(readset)>0:
                self.subpath.append((sofar1, readset))
            return
        
        # path not supported
        if len(readset) <1:
            return  
        
        # path not circular
        if start == "SINK":
            return
        
        depth1 = depth+ 1
        
        
        for dst in g.outgoing[start]:   
            # consider the read set in the latest 10 intervals
            if dst.startswith("E"):
                readset1 = readset & set(g.edges[dst]['reads'])
            else:
                readset1 = readset
            
            self.find_path(g, dst, end, sofar = sofar + [start], depth = depth1, readset = readset1)
            
def reconstruct_path_seq(graph, path):
    seq = ""
    for item in path:
        if item.startswith('A'):
            seq += graph.anchor[item]['seq']
        elif item.startswith("E"):
            seq += graph.edges[item]['seq']
        else:
            item += ""
    return seq

def find_furthest_node(node_candidate, subgraph):
    max_distance = 0
    node = ""
    for n in node_candidate:
        d = numpy.absolute(subgraph.anchor[node]['pos'] - subgraph.anchor['start_node']['pos'])
        if d > max_distance:
            node = n
    return node

class Get_Series_Parallel_Graph:
    
    def __init__(self, graph):
        self.initial_set = self.find_all_reads(graph)
        self.nodelist = self.series_parallel_graph_nodelist(graph)
        #print(self.nodelist)
        self.anchor, self.edges, self.outgoing, self.incoming = self.series_parallel_graph(self.nodelist, graph)

    def find_all_reads(self, graph):
        read_sets = set()
        edgelist = graph.edges.keys()
        for item in edgelist:
            readlist = graph.edges[item]['reads']
            for read in readlist:
                read_sets.add(read)
        return read_sets
    def find_furthest_node(self, node_candidate, subgraph, start_node):
        max_distance = 0
        node = ""
        for n in node_candidate:
            d = numpy.absolute(subgraph.anchor[n]['pos'] - subgraph.anchor[start_node]['pos'])
            if d > max_distance:
                node = n
        return node

    def series_parallel_graph_nodelist(self, subgraph): 
        #  need to consider the loop i.e.the last anchor to the first anchor

        start_node = sorted(subgraph.anchor.keys())[0]
        Nodelist = [start_node]

        edgelist = subgraph.outgoing[start_node]
        node_candidate = []
        for edge in edgelist:
            nodelist = subgraph.outgoing[edge]
            node_candidate += nodelist
            if nodelist[0] not in subgraph.anchor:
                    continue
        node_candidate = sorted(node_candidate)
        #node = node_candidate[-1] ## node should be selected by the largest distance instead of the sorting order
        node = self.find_furthest_node(node_candidate, subgraph, start_node)

        # find the furthurest anchor
        Nodelist.append(node) # append the furthest node
        


        while node != start_node:
            # print(node)

            edgelist = subgraph.outgoing[node]
            node_candidate = []
            for edge in edgelist:
                nodelist = subgraph.outgoing[edge]
                # exclude deadend
                if "SINK" in nodelist:
                    continue
                if nodelist[0] not in subgraph.anchor:
                    continue
                if nodelist[0] not in subgraph.outgoing:
                    continue
                node_candidate += nodelist

            node_candidate = sorted(node_candidate)
            node = self.find_furthest_node(node_candidate, subgraph, node)
            if node in set(Nodelist):
                Nodelist.append(node)
                break
            Nodelist.append(node) # append the furthest node  
                
        return Nodelist
    
    def series_parallel_graph(self, Nodelist, subgraph):
        Node_dict = {}
        Edge_dict = {}
        Outgoing_dict = {}
        Incoming_dict = {}
        for i, node in enumerate(Nodelist[:-1]):
            start_node = node
            end_node = Nodelist[i+1]
            Node_dict[start_node] = subgraph.anchor[start_node]
            # print(start_node, end_node)
            path = Find_all_Path_between_anchors(subgraph, start_node, end_node, self.initial_set)
            #print(start_node,end_node, len(path.subpath))
            index = 0
            for p, rs in path.subpath:
                edgename = 'E%05d.%04d' % (int(start_node[1:]), index)
                seq = reconstruct_path_seq(subgraph, p[1:-1])
                Edge_dict[edgename] = {}
                Edge_dict[edgename]['seq'] = seq
                Edge_dict[edgename]['src'] = start_node
                Edge_dict[edgename]['dst'] = end_node
                Edge_dict[edgename]['reads'] = list(rs)
                Edge_dict[edgename]['strain'] = list(set([item.split("_")[-1] for item in list(rs)]))

                Outgoing_dict[start_node] = Outgoing_dict.get(start_node, []) + [edgename]
                Outgoing_dict[edgename] = [end_node]

                Incoming_dict[end_node] = Incoming_dict.get(end_node, []) + [edgename]
                Incoming_dict[edgename] = [start_node]
                index += 1
        return Node_dict, Edge_dict, Outgoing_dict, Incoming_dict
    

class SubGraph:
    def __init__(self, graph, samplelist):
        self.anchor, self.edges, self.outgoing, self.incoming = self.reconstruct_sample_subgraph(graph, samplelist)

    def reconstruct_sample_subgraph(self, graph, samplelist):
        Anchor_dict = {}
        Edge_dict = {}
        Outgoing = {}
        Incoming = {}

        edgelist = list(graph.edges.keys())
        nodelist = []
        
        for edge in edgelist:
            if len(set(samplelist) & set(graph.edges[edge]['strain'])) > 0:
                Edge_dict[edge] = graph.edges[edge]
                src = graph.incoming[edge][0]
                dst = graph.outgoing[edge][0]
                
                Edge_dict[edge]['src'] = src
                Edge_dict[edge]['dst'] = dst
                
                Incoming[edge] = graph.incoming[edge]
                Incoming[dst] = Incoming.get(dst, []) + [edge]
                
                Outgoing[edge] = graph.outgoing[edge]
                Outgoing[src] = Outgoing.get(src, []) + [edge]
                
                nodelist += graph.incoming[edge]
                nodelist += graph.outgoing[edge]

        nodelist = list(set(nodelist))
        for node in nodelist:
            if node.startswith("A"):
                Anchor_dict[node] = graph.anchor[node]      
        return Anchor_dict, Edge_dict, Outgoing, Incoming
    
def validate_sp_graph(series_parallelgraph):
    nodelist = series_parallelgraph.anchor.keys()
    for node in nodelist:
        edgelist = series_parallelgraph.outgoing[node]
        outgoing_nodelist = []
        for edge in edgelist:
            outgoing_nodelist.append(series_parallelgraph.edges[edge]['dst'])
            outgoing_nodelist.append(series_parallelgraph.outgoing[edge][0])
        outgoing_nodelist = set(outgoing_nodelist)
        assert len(outgoing_nodelist) == 1
    print("PASS")

def processCigar(cigar):
    """Helper Function, may not be used directly, expand Cigar string
    
    Parameters:
        cigar: <str> - compressed cigar
    """
    out = ''
    N = 0
    for symbol in cigar:
        if symbol in '0123456789':
            N = 10*N + int(symbol)
        else:
            #if (symbol != 'D'):
            if (N == 0):
                out += symbol
            else:
                out += N*symbol
            N = 0
    return out

def combineCigar(cigar):
    """Helper Function, may not be used directly, compress Cigar string
    
    Parameters:
        cigar: <str> - expanded cigar
    """
    cigar = cigar +'$'
    out = ''
    N = 0
    start = 0
    for i in range(1,len(cigar)):
        if cigar[i-1] != cigar[i]:
            out += str(i-start) + cigar[i-1]
            start = i
    return out  

def get_variant_position(cigar):
    ref_pos = []
    alt_pos = []
    var_type = []
    alt_i = 0
    ref_i = 0
    for i, s in enumerate(cigar):
        if s == 'I':
            if ref_i > 0:
                ref_pos.append(ref_i-1)
            else:
                ref_pos.append(ref_i)    
            alt_pos.append(alt_i)
            var_type.append("I")
            
            alt_i += 1
            
        if s == 'D':
            ref_pos.append(ref_i)
            if alt_i > 0:
                alt_pos.append(alt_i-1)
            else:
                alt_pos.append(alt_i)
            var_type.append("D")
            ref_i += 1
            
        if s == 'X':
            ref_pos.append(ref_i)
            alt_pos.append(alt_i)
            alt_i += 1
            ref_i += 1
            var_type.append("X")
            
        if s == '=':
            alt_i += 1
            ref_i += 1

    return ref_pos, alt_pos, var_type

def findBedge(Graph, src, dst, refstrain, k):
    paths = AGG.Find_all_Path_between_anchors(Graph, src, dst, {refstrain})
    subpaths = paths.subpath

    if len(subpaths) < 1:
        return ""
    for p, strain in subpaths:
        seq = AGG.reconstruct_path_seq(spgraph, path = p)
    return seq[k:-k]

def find_all_reads(graph):
    read_sets = set()
    edgelist = graph.edges.keys()
    for item in edgelist:
        readlist = graph.edges[item]['reads']
        for read in readlist:
            read_sets.add(read)
    return read_sets

def get_snps(subgraph, mpath, k, ref):
    SNPs = {} # Var[refpos]['edgename'] = ['A']
    for edge in mpath:
        if edge.startswith('E'):
            cigar = subgraph.edges[edge].get('variants', "")
            if cigar == "":
                continue
            src = subgraph.incoming[edge][0]
            if src != "SOURCE":
                refstart = int(subgraph.anchor[src]["pos"]) + k
            else:
                raise "SOURCE or SINK node"

            # find reference seq
            dst = subgraph.outgoing[edge][0]
            ref_edge_list, ref_seq = findBedge(spgraph, src, dst, "NC_012920", k)
            ref_edge = ref_edge_list[0]
            expanded_cigar = AGG.processCigar(cigar)
            refpos, altpos, var_type = get_variant_position(AGG.processCigar(cigar))
            alt_seq = subgraph.edges[edge]['seq']
            # only record SNPs
            for i, vart in enumerate(var_type):
                if vart == "X":
                    rp = refpos[i]
                    ap = altpos[i]
                    pos = refstart + rp
                    SNPs[pos] = SNPs.get(pos, {})
                    SNPs[pos][edge] = SNPs[pos].get(edge, {})                
                    SNPs[pos][edge]['base'] = SNPs[pos][edge].get('base', "") + alt_seq[ap]
                    SNPs[pos][edge]['altoffset'] = SNPs[pos][edge].get('altoffset', []) + [ap]
                    SNPs[pos][ref_edge] = SNPs[pos].get(ref_edge, {})
                    SNPs[pos][ref_edge]['base'] = ref_seq[refpos[i]]
                    SNPs[pos]['refbase'] = SNPs[pos][ref_edge]['base']
    return SNPs

def get_indels(subgraph, mpath, k, ref):
    Indels = {} # Var[refpos]['edgename'] = ['A']
    for edge in mpath:
        if edge.startswith('E'):
            cigar = subgraph.edges[edge].get('variants', "")
            if cigar == "":
                continue
            src = subgraph.incoming[edge][0]
            if src != "SOURCE":
                refstart = int(subgraph.anchor[src]["pos"]) + k
            else:
                raise "SOURCE or SINK node"

            # find reference seq
            dst = subgraph.outgoing[edge][0]
            ref_edge_list, ref_seq = findBedge(spgraph, src, dst, "NC_012920", k)
            ref_edge = ref_edge_list[0]
            expanded_cigar = AGG.processCigar(cigar)
            refpos, altpos, var_type = get_variant_position(AGG.processCigar(cigar))
            alt_seq = subgraph.edges[edge]['seq']
            # only record SNPs
            for i, vart in enumerate(var_type):
                if vart == "I":
                    rp = refpos[i]
                    ap = altpos[i]
                    pos = refstart + rp
                    Indels[pos] = Indels.get(pos, {})
                    Indels[pos][edge] = Indels[pos].get(edge, {})                
                    Indels[pos][edge]['base'] = Indels[pos][edge].get('base', "") + alt_seq[ap:ap+2]
                    Indels[pos][edge]['info'] = "INS"
                    Indels[pos][edge]['altoffset'] = Indels[pos][edge].get('altoffset', []) + [ap]
                    Indels[pos][ref_edge] = Indels[pos].get(ref_edge, {})
                    Indels[pos][ref_edge]['base'] = Indels[pos][ref_edge].get("base", "") + ref_seq[rp-1]
                    Indels[pos]['refbase'] = Indels[pos][ref_edge]['base']

                if vart == "D":
                    rp = refpos[i]
                    ap = altpos[i]
                    pos = refstart + rp
                    Indels[pos] = Indels.get(pos, {})
                    Indels[pos][edge] = Indels[pos].get(edge, {})                
                    Indels[pos][edge]['base'] = Indels[pos][edge].get('base', "") + alt_seq[ap]
                    Indels[pos][edge]['info'] = "DEL"
                    Indels[pos][edge]['altoffset'] = Indels[pos][edge].get('altoffset', []) + [ap-1]
                    Indels[pos][ref_edge] = Indels[pos].get(ref_edge, {})
                    Indels[pos][ref_edge]['base'] = Indels[pos][ref_edge].get("base", "") + ref_seq[rp-1:rp+1]
                    Indels[pos]['refbase'] = Indels[pos][ref_edge]['base']
    return Indels


def get_VCF_file(Var, samplelist):
    Mat = {}

    for pos, D in Var.items():
        Mat[pos] = {}
        for edge, B in D.items():
            if edge == "refbase":
                continue
            for strain in samplelist:
                Mat[pos][strain] = Mat[pos].get(strain, set()) | set([B['base']])
    #print(Mat)
    # write VCF
    VCF = {}
    for pos, D in Mat.items():
        refbase = Var[pos]['refbase']
        if list(D.values())[0] == set([refbase]):
            continue    
        VCF[pos]= VCF.get(pos, {})
        for strain, base in D.items():
            VCF[pos]["CHROM"] = "chrM"
            VCF[pos]["POS"] = pos
            VCF[pos]['ID'] = "."
            VCF[pos]['REF'] = refbase
            #print(base - set([refbase]))
            VCF[pos]['ALT'] = ",".join(base - set([refbase]))
            VCF[pos]['QUAL'] = "."
            VCF[pos]['FILTER'] = "."
            VCF[pos]['INFO'] = "."
            VCF[pos]['FORMAT'] = "GT"
            #print(base)
            # allelelist = []
            # if refbase in base:
            #     allelelist += [0]

            # for i,s in enumerate(base - set([refbase])):
            #     allelelist += [i+1]
            # VCF[pos][strain] = "/".join([str(item) for item in allelelist])
            VCF[pos][strain] = "/".join([str(item) for item in base])
    return VCF


def loadFasta(filename):
    """ Parses a classically formatted and possibly 
        compressed FASTA file into a list of headers 
        and fragment sequences for each sequence contained.
        The resulting sequences are 0-indexed! """
    if (filename.endswith(".gz")):
        fp = gzip.open(filename, 'rb')
    else:
        fp = open(filename, 'rb')
    # split at headers
    data = fp.read().decode().split('>')
    fp.close()
    # ignore whatever appears before the 1st header
    data.pop(0)     
    headers = []
    sequences = []
    for sequence in data:
        lines = sequence.split('\n')
        headers.append(lines.pop(0))
        sequences.append(''.join(lines))
    return (headers, sequences)