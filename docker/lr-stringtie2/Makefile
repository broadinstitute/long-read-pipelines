VERSION = 2.2.1-test # This should match the Stringtie version number

TAG1 = us.gcr.io/broad-dsp-lrma/lr-stringtie2:$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/lr-stringtie2:latest

all: | build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

build_no_cache:
	docker build --no-cache -t $(TAG1) -t $(TAG2) .
  
push:
	docker push $(TAG1)
	docker push $(TAG2)
