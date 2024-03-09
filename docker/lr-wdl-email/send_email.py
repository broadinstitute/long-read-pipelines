import argparse
import base64
import os
from typing import List
import pandas as pd

from sendgrid import SendGridAPIClient
from sendgrid.helpers.mail import (Attachment, FileContent, FileName, FileType, Disposition, ContentId)
from sendgrid.helpers.mail import From as SGFrom, To as SGTo, Subject as SGSubject
from sendgrid.helpers.mail import Mail, PlainTextContent, HtmlContent


def send_notification(notification_sender_name: str,
                      notification_receiver_names: List[str],
                      notification_receiver_emails: List[str],
                      email_subject: str,
                      email_body: str,
                      html_body: str = None) -> None:
    """
    Sending notification email to (potentially) multiple recipients.

    Note that this assumes two environment variables are set appropriately: "SENDGRID_API_KEY" & "SENDER_EMAIL".

    Provide html_body at your own risk.

    Shameless copy from
    https://github.com/sendgrid/sendgrid-python/blob/main/examples/helpers/mail_example.py#L9
    :return:
    """
    if len(notification_receiver_emails) != len(notification_receiver_names):
        raise ValueError("Different number of recipients and recipients' emails")

    sg = SendGridAPIClient(api_key=os.environ.get('SENDGRID_API_KEY'))

    email_core = _construct_sendgrid_mail_core(notification_sender_name, email_subject, email_body, html_body)

    # send
    failed_responses = list()
    from copy import deepcopy
    for i in range(len(notification_receiver_emails)):
        message = deepcopy(email_core)
        message.add_to(SGTo(notification_receiver_emails[i], notification_receiver_names[i]))
        response = sg.client.mail.send.post(request_body=message.get())
        if 202 != response.status_code:
            failed_responses.append(i)
    if 0 < len(failed_responses):
        failures = '  \n'.join([notification_receiver_names[i]+':'+notification_receiver_emails[i]
                                for i in failed_responses])
        logger.warning(f"Failed to send message to some receivers: \n  {failures}")


def send_notification_with_attachments(notification_sender_name: str,
                                       notification_receiver_names: List[str],
                                       notification_receiver_emails: List[str],
                                       email_subject: str,
                                       email_body: str,
                                       html_body: str = None,
                                       txt_names_and_files: list = None,
                                       tsv_names_and_files: list = None,
                                       pdf_names_and_files: list = None
                                       ) -> None:
    """
    Sending notification email to (potentially) multiple recipients.
    Note that this assumes two environment variables are set appropriately, "SENDGRID_API_KEY" & "SENDER_EMAIL".

    Provide html_body at your own risk.

    Shameless copy from
    https://github.com/sendgrid/sendgrid-python/blob/main/examples/helpers/mail_example.py#L9
    :return:
    """
    if len(notification_receiver_emails) != len(notification_receiver_names):
        raise ValueError("Different number of recipients and recipients' emails")

    sg = SendGridAPIClient(api_key=os.environ.get('SENDGRID_API_KEY'))

    email_core = _construct_sendgrid_mail_core(notification_sender_name, email_subject, email_body, html_body)

    # attach
    attachments = _attach_files_to_mail(txt_names_and_files, tsv_names_and_files, pdf_names_and_files)
    for a in attachments:
        email_core.add_attachment(a)

    # send
    failed_responses = list()
    from copy import deepcopy
    for i in range(len(notification_receiver_emails)):
        message = deepcopy(email_core)
        message.add_to(SGTo(notification_receiver_emails[i], notification_receiver_names[i]))
        response = sg.client.mail.send.post(request_body=message.get())
        if 202 != response.status_code:
            failed_responses.append(i)
    if 0 < len(failed_responses):
        failures = '  \n'.join([notification_receiver_names[i]+':'+notification_receiver_emails[i]
                                for i in failed_responses])
        logger.warning(f"Failed to send message to some receivers: \n  {failures}")


def _construct_sendgrid_mail_core(notification_sender_name: str,
                                  email_subject: str,
                                  email_body: str,
                                  html_body: str = None) -> Mail:

    """
    Construct core content of notification email to (potentially) multiple recipients.

    Note that this assumes two environment variables are set appropriately, "SENDGRID_API_KEY" & "SENDER_EMAIL".

    The returned Mail object DOES NOT specify recipients, caller should customize that.

    Provide html_body at your own risk.

    Shameless copy from
    https://github.com/sendgrid/sendgrid-python/blob/main/examples/helpers/mail_example.py#L9
    :return:
    """

    # arg validation
    assert "SENDGRID_API_KEY" in os.environ, \
        'environment variable SENDGRID_API_KEY is needed.'
    assert "SENDER_EMAIL" in os.environ, \
        'environment variable SENDER_EMAIL is needed.'

    # construct core
    notification_sender_email = os.environ.get('SENDER_EMAIL')
    email_core = Mail(from_email=SGFrom(notification_sender_email, notification_sender_name),
                      subject=SGSubject(email_subject),
                      plain_text_content=PlainTextContent(email_body),
                      html_content=HtmlContent(html_body) if html_body else None,
                      is_multiple=True)  # recipients won't see each other
    return email_core


def _attach_files_to_mail(txt_names_and_files: list = None,
                          tsv_names_and_files: list = None,
                          pdf_names_and_files: list = None) -> list:
    """
    Return a list of SendGrid Attachments for attaching to the email core.

    :param txt_names_and_files: list of tuple2 (file name, file path)
    :param tsv_names_and_files: list of tuple2 (file name, file path)
    :param pdf_names_and_files: list of tuple2 (file name, file path)
    :return: a list of SendGrid Attachments for attaching to the email core.
    """

    has_something_to_attach = False
    for a in [txt_names_and_files, tsv_names_and_files, pdf_names_and_files]:
        if a is not None and 0 != len(a):
            has_something_to_attach = True
            break

    assert has_something_to_attach, "No valid inputs for building attachments"

    attachments = list()

    # txt
    if txt_names_and_files is not None:
        for attachment_txt_name, attachment_txt_path in txt_names_and_files:
            with open(attachment_txt_path) as inf:
                contents = [line.strip() for line in inf.readlines()]
            base64_txt = \
                base64.b64encode(('\n'.join(contents)).encode('utf-8')).decode('utf-8')
            txt_attachment = Attachment(
                FileContent(base64_txt),
                FileName(attachment_txt_name),
                FileType('text/plain'),
                Disposition('attachment'),
                ContentId(attachment_txt_name)
            )
            attachments.append(txt_attachment)
    # tsv
    if tsv_names_and_files is not None:
        for attachment_tsv_name, attachment_tsv_path in tsv_names_and_files:
            temp_dataframe = pd.read_csv(attachment_tsv_path, header='infer', sep='\t')
            base64_csv = \
                base64.b64encode(temp_dataframe.to_csv(header=True, index=False, sep='\t').encode()).decode()
            tsv_attachment = Attachment(
                FileContent(base64_csv),
                FileName(attachment_tsv_name),
                FileType('text/csv'),
                Disposition('attachment'),
                ContentId('dataframe')
            )
            attachments.append(tsv_attachment)
    # pdf
    if pdf_names_and_files is not None:
        for attachment_pdf_name, attachment_pdf_path in pdf_names_and_files:
            with open(attachment_pdf_path, 'rb') as f:
                data = f.read()
            pdf_content = base64.b64encode(data).decode()

            pdf_attachment = Attachment(
                FileContent(pdf_content),
                FileName(attachment_pdf_name),
                FileType('application/pdf'),
                Disposition('attachment')
            )
            attachments.append(pdf_attachment)

    return attachments


def main():
    parser = argparse.ArgumentParser(description='Send email via SendGrid',
                                     prog='send_email')

    parser.add_argument('--sendgrid_api_key', type=str,
                        help="JSON file holding the secret API key registered at SendGrid")
    parser.add_argument('--sender_name', type=str,
                        help="Name of the sender, i.e. identify yourself")
    parser.add_argument('--sender_email', type=str,
                        help="Email address registered at SendGrid for sending out the email")

    parser.add_argument('--notification_receiver_names',
                        type=str,
                        help="Read names of duplicate records")
    parser.add_argument('--notification_receiver_emails',
                        type=str,
                        help="Read names of duplicate records")

    parser.add_argument('--email_subject', type=str,
                        help="The subject/title/topic of the email")
    parser.add_argument('--email_body', type=str,
                        help="The plain-text contents of the email")

    # we hide this option in cli because that's for later
    # parser.add_argument('--html_body', type=str, help="if you")

    parser.add_argument('--txt_names_and_files', type=str,
                        help="2-col TSV holding the (desired name, path) of the txt file to be attached")
    parser.add_argument('--tsv_names_and_files', type=str,
                        help="2-col TSV holding the (desired name, path) of the tsv file to be attached")
    parser.add_argument('--pdf_names_and_files', type=str,
                        help="2-col TSV holding the (desired name, path) of the pdf file to be attached")

    args = parser.parse_args()

    with open(args.sendgrid_api_key, 'r') as inf:
        os.environ['SENDGRID_API_KEY'] = inf.readlines()[0].strip()
    os.environ['SENDER_EMAIL'] = args.sender_email

    with open(args.notification_receiver_names) as inf:
        receivers = [n.strip() for n in inf.readlines()]
    with open(args.notification_receiver_emails) as inf:
        receiver_emails = [n.strip() for n in inf.readlines()]

    def load_from_headerless_tsv(ff: str) -> list or None:
        if ff is None:
            return None
        else:
            a = pd.read_csv(ff, header=None, sep='\t')
            return list(zip(a.iloc[:, 0], a.iloc[:, 1]))

    send_notification_with_attachments(notification_sender_name=args.sender_name,
                                       notification_receiver_names=receivers,
                                       notification_receiver_emails=receiver_emails,

                                       email_subject=args.email_subject,
                                       email_body=args.email_body,
                                       html_body=None,

                                       txt_names_and_files=load_from_headerless_tsv(args.txt_names_and_files),
                                       tsv_names_and_files=load_from_headerless_tsv(args.tsv_names_and_files),
                                       pdf_names_and_files=load_from_headerless_tsv(args.pdf_names_and_files))


if __name__ == "__main__":
    main()
