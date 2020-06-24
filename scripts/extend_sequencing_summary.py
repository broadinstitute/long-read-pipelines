import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-s', dest='sequencing_summary', type=str,
                    help='Sequencing Summary')
parser.add_argument('-r', dest='read_ids', type=str,
                    help='Read IDs')
args = parser.parse_args()

read_ids = {}

for line in open(args.read_ids, 'r'):
    read_id = line.strip()
    original_read_id = read_id.split('_')[0]

    if original_read_id not in read_ids:
        read_ids[original_read_id] = []

    read_ids[original_read_id].append(read_id)

read_id_idx = None

for line in open(args.sequencing_summary, 'r'):
    tokens = map(lambda s: s.strip(), line.split('\t'))
    if read_id_idx is None:
        read_id_idx = tokens.index('read_id')
        print("\t".join(tokens))
        continue

    #For a given line in the seq summary. print it n number of times for the size of the read_ids[read_id] array while swapping out the read_id
    for new_read_id in read_ids[tokens[read_id_idx]]:
        tokens[read_id_idx] = new_read_id
        print("\t".join(tokens))
