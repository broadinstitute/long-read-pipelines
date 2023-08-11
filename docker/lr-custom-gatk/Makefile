VERSION = 4.4.0.0-samtools1.18
TAG1 = us.gcr.io/broad-dsp-lrma/lr-custom-gatk:$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/lr-custom-gatk:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
