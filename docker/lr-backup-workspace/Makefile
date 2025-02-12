VERSION = 0.0.1
TAG1 = us.gcr.io/broad-dsp-lrma/lr-backup-workspace:$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/lr-backup-workspace:latest

all: build
 
build:
	docker build --platform linux/amd64	-t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
