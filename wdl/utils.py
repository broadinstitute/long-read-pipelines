from subprocess import check_output

def list_jobs():
    # get the last N worflows
    jobs = check_output("./list_jobs.py")
    lines = jobs.decode('UTF-8').strip().split('\n')
    records = []
    for line in lines:
        if "broadinstitute.org" not in line:
            continue
        fields = line.split()
        date_string = fields[0]
        workflow_id = fields[2]
        name = fields[3]
        status = fields[4]
        
        records.append({
                'date': date_string,
                'workflow_id': workflow_id,
                'name': name,
                'status': status,
        })
    
    return records

def list_workflow_ids():
    return  [r["workflow_id"] for r in list_jobs()]

def parse_workflow_ids(workflow_ids):
    if len(workflow_ids) == 1:
        try:
            num_lines = abs(int(workflow_ids[0]))
        except ValueError:
            pass
        else:
            # get the last N worflows
            workflow_ids = list_workflow_ids()[-1*num_lines:]

    return workflow_ids
