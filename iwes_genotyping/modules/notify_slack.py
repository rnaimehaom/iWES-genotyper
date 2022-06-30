import time
from slack import WebClient
from slack.errors import SlackApiError
import sys
from dotenv import load_dotenv
import os
load_dotenv()

#token = os.environ.get("api-token")

program_name = sys.argv[0]
cassowary_path = sys.argv[1]
submit_name = sys.argv[2]


slack_message = '{0} finished and is in: \n{1}'.format(submit_name, cassowary_path)

client = WebClient(token=os.environ.get('SLACK_BOT_TOKEN2'))

user_list = ['UA7CAHQH2',
             'U029UCMM5E1']

for user_i in user_list:
    try:
        response = client.chat_postMessage(
            channel=user_i,
            text=slack_message)
        assert response["message"]["text"] == slack_message
        print('notified {0}'.format(user_i))
    except SlackApiError as e:
        # You will get a SlackApiError if "ok" is False
        assert e.response["ok"] is False
        assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
        print(f"Got an error: {e.response['error']}")
    time.sleep(5)

