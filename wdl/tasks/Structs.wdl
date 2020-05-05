version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    String? disk_type
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}
