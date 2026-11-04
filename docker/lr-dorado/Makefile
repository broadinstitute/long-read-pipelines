IMAGE_NAME = lr-dorado
# The bundled Dorado binary is 1.4.0; this image rev (1.4.1) adds the pod5 CLI used by the
# fast5 -> pod5 conversion step, so the tag no longer exactly equals the Dorado version.
VERSION = 1.4.1

TAG1 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME):$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME):latest

all: | build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

build_no_cache:
	docker build --no-cache -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
