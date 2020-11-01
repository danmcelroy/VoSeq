import logging
from smtplib import SMTPException

from django.conf import settings
from django.contrib.auth.models import User
from django.core.mail import send_mail
from django.http import HttpRequest
from django.urls import reverse

from voseq.celery import app

log = logging.getLogger(__name__)


@app.task
def log_email_error(
    request: HttpRequest, exc: str, traceback: str, task_id: str
) -> None:
    log.error(
        f"log_email_error\n--\n\nrequest {request} \n\nexc {exc}"
        f"\n\ntraceback {traceback}"
    )


@app.task
def notify_user(job, user_id) -> None:
    """Send an email notification to user

    If the gui_user is specified, we will send the notification to the person
    that is doing actions via the GUI. Otherwise, we will notify the user that
    created the ContactJob.
    """
    gui_user = User.objects.get(id=user_id)
    user = gui_user if gui_user else job.user
    log.debug("notify_user " + str(job))

    subject = "Person Job Completed - " + str(job.name)
    relative_url = reverse('create_dataset.results', args=(job.id,))
    result_url = "http://app.hiplead.com" + relative_url
    content = "Your contact job {0} has successfully completed. " \
              "Please verify and download the results from: " \
              "{1}".format(job.name, result_url)
    from_email = 'noreply@hiplead.com'

    if user and user.email:
        to_emails = [user.email] + [email for name, email in settings.ADMINS]
        try:
            send_mail(subject, content, from_email, to_emails)
        except SMTPException:
            log.exception("Failed to notify_list_users for job " + str(job))
        else:
            log.debug("sent list change status email to " + str(to_emails))
    else:
        log.debug('Cannot send notification email. '
                  'No user / email assigned to job ' + str(job))