'''app.py
Building the slack app
July 1, 2021
author: Johanna Harding
helpful sources: https://github.com/slackapi/bolt-python/tree/main/examples'''

import os
# Use the package we installed
from flask import Flask, redirect, url_for
# from .flask import app
from slack_bolt import App
# import recreation #will need recreation at some point, but leaving it out for now



app = Flask(__name__)


os.environ['SLACK_BOT_TOKEN'] = 'xoxb-2045402657443-2250384640595-I9MhN3qoyPcVO92cgloUdokE'
os.environ['SLACK_SIGNING_SECRET'] = ''

sl_app = App(
  token = os.getenv('SLACK_BOT_TOKEN'),
  secret = os.getenv('SLACK_SIGNING_SECRET')
)


@app.route('app_home_opened')
def app_home_tab(client, event):
  '''want to build the home tab of the app'''
  client.views_publish(
    user_id=event["user"],
    view={
      "type": "home",
      "callback_id": "home_view",
       "blocks": [
          {
            "type": "section",
            "text": {
              "type": "mrkdwn",
              "text": "*Welcome to Observation Automation*"
            }
          },
          {
            "type": "divider"
          },
          {
            "type": "section",
            "text": {
              "type": "mrkdwn",
              "text": "This button won't do much for now but you can set up a listener for it using the `actions()` method and passing its unique `action_id`. See an example in the `examples` folder within your Bolt app."
            }
          },
          {
            "type": "actions",
            "elements": [
              {
                "type": "button",
                "text": {
                  "type": "plain_text",
                  "text": "Submit an Observation",
                  "action_id": "submit_observation"
                }
              }
            ]
          },
          {
            "type": "actions",
            "elements":[
              {
                "type": "button",
                "text": {
                  "type": "plain_text",
                  "text": "Check Observation Status"
                }
              }
            ]
          }
        ]
      })


@app.route('/')#, methods['POST'])
# @login_required
def submit_observation_action(ack):
  ack()
  redirect('/submit_observation')
  # if request.method == 'POST':
  #   print('here')

@app.route('/submit_observation')
def new_submission():
  '''should give options for a new submission'''
  {
	"type": "",
	"blocks": [
    {
			"dispatch_action": true,
			"type": "input",
			"element": {
				"type": "plain_text_input",
				"action_id": "plain_text_input-action"
			},
			"label": {
				"type": "plain_text",
				"text": "Label",
				"emoji": true
			}
    },
    {
			"dispatch_action": true,
			"type": "input",
			"element": {
				"type": "plain_text_input",
				"action_id": "plain_text_input-action"
			},
			"label": {
				"type": "plain_text",
				"text": "Label",
				"emoji": true
      }
    }
  ]
  }



# app = App(
#     token=os.environ.get("SLACK_BOT_TOKEN"),
#     signing_secret=os.environ.get("SLACK_SIGNING_SECRET")
# )

# # Add functionality here
# # @app.event("app_home_opened") etc
# @app.event("app_home_opened")
# def update_home_tab(client, event, logger):
#   try:
#     # views.publish is the method that your app uses to push a view to the Home tab
#     client.views_publish(
#       # the user that opened your app's app home
#       user_id=event["user"],
#       # the view object that appears in the app home
#       view={
#         "type": "home",
#         "callback_id": "home_view",

#         # body of the view
#         "blocks": [
#           {
#             "type": "section",
#             "text": {
#               "type": "mrkdwn",
#               "text": "*Welcome to your _App's Home_* :tada:"
#             }
#           },
#           {
#             "type": "divider"
#           },
#           {
#             "type": "section",
#             "text": {
#               "type": "mrkdwn",
#               "text": "This button won't do much for now but you can set up a listener for it using the `actions()` method and passing its unique `action_id`. See an example in the `examples` folder within your Bolt app."
#             }
#           },
#           {
#             "type": "actions",
#             "elements": [
#               {
#                 "type": "button",
#                 "text": {
#                   "type": "plain_text",
#                   "text": "Click me!"
#                 }
#               }
#             ]
#           }
#         ]
#       }
#     )
  
#   except Exception as e:
#     logger.error(f"Error publishing home tab: {e}")


# # Start your app
# if __name__ == "__main__":
#     app.start(port=int(os.environ.get("PORT", 3000)))
