#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This file is part of CMInject
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# If you use this program for scientific work, you should correctly reference it; see LICENSE file for details.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not, see
# <http://www.gnu.org/licenses/>.

from abc import ABC, abstractmethod
from collections import defaultdict
from enum import Enum
from typing import Any, Union, List, Iterable, Callable, Dict, Set

from cminject.utils.software_structure import Singleton


__all__ = ['ConfigKey', 'ConfigSubscriber', 'GlobalConfig']


class ConfigKey(Enum):
    NUMBER_OF_DIMENSIONS = 1
    TIME_STEP = 2


class ConfigSubscriber(ABC):
    @abstractmethod
    def config_change(self, key: ConfigKey, value: Any):
        pass


@Singleton
class GlobalConfig:
    def __init__(self):
        self._subscribers: Dict[ConfigKey, Set[ConfigSubscriber]] = defaultdict(set)
        self._config: Dict[ConfigKey, Any] = {k: None for k in ConfigKey}

    def set(self, key: ConfigKey, value: Any):
        self._config[key] = value
        for subscriber in self._subscribers.get(key, []):
            subscriber.config_change(key, self._config[key])

    def get(self, key: ConfigKey) -> Any:
        return self._config[key]

    @staticmethod
    def _call_for_all_pairs(subscriber: Union[ConfigSubscriber, Iterable[ConfigSubscriber]],
                            key: Union[ConfigKey, List[ConfigKey]],
                            fn: Callable[[ConfigSubscriber, ConfigKey], Any]):
        if not isinstance(subscriber, Iterable):
            subscriber = [subscriber]
        if not isinstance(key, list):
            key = [key]

        for key_ in key:
            for subscriber_ in subscriber:
                fn(subscriber_, key_)

    def _subscribe_one(self, subscriber: ConfigSubscriber, key: ConfigKey):
        self._subscribers[key].add(subscriber)
        # notify once at the moment of subscription, if value for this key has already been set
        if self._config[key] is not None:
            subscriber.config_change(key, self._config[key])

    def subscribe(self, subscriber: Union[ConfigSubscriber, Iterable[ConfigSubscriber]],
                  key: Union[ConfigKey, List[ConfigKey]]):
        return self._call_for_all_pairs(subscriber, key, self._subscribe_one)

    def _unsubscribe_one(self, subscriber: ConfigSubscriber, key: ConfigKey):
        self._subscribers[key].remove(subscriber)

    def unsubscribe(self, subscriber: Union[ConfigSubscriber, Iterable[ConfigSubscriber]],
                    key: Union[ConfigKey, List[ConfigKey]]):
        return self._call_for_all_pairs(subscriber, key, self._unsubscribe_one)

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
