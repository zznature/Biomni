import os

from langchain_anthropic import ChatAnthropic
from langchain_openai import AzureChatOpenAI, ChatOpenAI


def get_llm(
    model="claude-3-5-sonnet-20241022",
    temperature=0.7,
    stop_sequences=None,
    azure=False,
):
    if azure:
        API_VERSION = "2024-12-01-preview"
        return AzureChatOpenAI(
            openai_api_key=os.getenv("OPENAI_API_KEY"),
            azure_endpoint=os.getenv("OPENAI_ENDPOINT"),
            azure_deployment=model,
            openai_api_version=API_VERSION,
        )

    if model[:7] == "claude-":
        source = "Anthropic"
    elif model[:4] == "gpt-":
        source = "OpenAI"
    if source not in ["OpenAI", "Anthropic"]:
        raise ValueError("Invalid source")
    if source == "OpenAI":
        return ChatOpenAI(model=model, temperature=temperature, stop_sequences=stop_sequences)
    elif source == "Anthropic":
        return ChatAnthropic(
            model=model,
            temperature=temperature,
            max_tokens=8192,
            stop_sequences=stop_sequences,
        )
