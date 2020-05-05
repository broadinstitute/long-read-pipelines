IMAGE_NAME = lr-jupyter
VERSION = 0.0.3

INTERACTIVE_SUFFIX = _interactive

TAG1 = quay.io/broad-long-read-pipelines/$(IMAGE_NAME):$(VERSION)
TAG2 = quay.io/broad-long-read-pipelines/$(IMAGE_NAME):latest
TAG3 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME):$(VERSION)
TAG4 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME):latest

TAG5 = quay.io/broad-long-read-pipelines/$(IMAGE_NAME)$(INTERACTIVE_SUFFIX):$(VERSION)
TAG6 = quay.io/broad-long-read-pipelines/$(IMAGE_NAME)$(INTERACTIVE_SUFFIX):latest
TAG7 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME)$(INTERACTIVE_SUFFIX):$(VERSION)
TAG8 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME)$(INTERACTIVE_SUFFIX):latest

all: | build push

build:
	docker build -t $(TAG1) -t $(TAG2) -t $(TAG3) -t $(TAG4) .
	cat Dockerfile | sed 's#CMD \["start-notebook.sh"\]#CMD \["/bin/bash"\]#g' > Dockerfile.interactive
	docker build - -t $(TAG5) -t $(TAG6) -t $(TAG7) -t $(TAG8) . < Dockerfile.interactive

build_no_cache:
	docker build --no-cache -t $(TAG1) -t $(TAG2) -t $(TAG3) -t $(TAG4) .
	cat Dockerfile | sed 's#CMD \["start-notebook.sh"\]#CMD \["/bin/bash"\]#g' > Dockerfile.interactive
	docker build - -t $(TAG5) -t $(TAG6) -t $(TAG7) -t $(TAG8) . < Dockerfile.interactive

push:
	docker push $(TAG1)
	docker push $(TAG2)
	docker push $(TAG3)
	docker push $(TAG4)
	docker push $(TAG5)
	docker push $(TAG6)
	docker push $(TAG7)
	docker push $(TAG8)

